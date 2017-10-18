clc;
clear;
close all;

%%% Error Metric Optimization

% This script computes the error metric minimum for a given number of
% robots, N, and prescribed robot radius, delta.

% To use this script, the user must enter the parameters noted by the USER
% INPUT sections below. These inputs include the swarm parameters, the
% computational parameters, and the target distribution.

% USER INPUT - Looping Parameter:
nLoop = 50;         % number of times to run the optimization

% Check for Existing Optimization Results:
errorsExist = exist('errors.mat','file');
if errorsExist == 0
    errors = zeros(1,nLoop);
    ki = 1;
    kf = nLoop;
    save('errors.mat','errors');
else
    load('errors.mat','errors');
    ki = length(errors) + 1;
    kf = length(errors) + nLoop;
    errors = [errors,zeros(1,nLoop)];
end

% Main Optimization Loop:
for k = ki:kf
    try
                
        disp(['Optimization k = ',num2str(k),' started.']); % start message
        
        % USER INPUT - Swarm and Computational Parameters:
        N = 200;               	% number of robots
        delta = 0.2;            % effective robot radius
        xmin = 0;               % left boundary
        xmax = 4.8;             % right boundary
        ymin = 0;               % bottom boundary
        ymax = 7;               % top boundary
        nx = 50;                % number of grid points, x
        ny = 50;                % number of grid points, y
        
        % Computed Parameters:
        xmid = (xmin+xmax)/2;   % center, x
        ymid = (ymin+ymax)/2;   % center, y
        dx = (xmax-xmin)/nx;    % grid spacing, x
        dy = (ymax-ymin)/ny;    % grid spacing, y
        X = xmin:dx:xmax;       % x grid
        Y = ymin:dy:ymax;       % y grid
        [X,Y] = meshgrid(X,Y);  % grid
        
        % USER INPUT - Target Distribution:
        % Note: The sample provided below is a ring target distribution.
        outerRadius = 2.06;  	% outer radius of ring
        innerRadius = 1.14;  	% inner radius of ring
        rho1 = 1;               % target distribution value, r<innerRadius
        rho2 = 36;              % target distribution value, innerRadius<r<outerRadius
        rho3 = rho1;            % target distribution value, r>outerRadius
        rho = (rho2-rho1)*double(((X-xmid).^2 ...
            + (Y-ymid).^2 < outerRadius^2) & ...
            ((X-xmid).^2 + (Y-ymid).^2 > innerRadius^2))...
            + rho1;             % target distribution definition
        rhoNorm = pi*(innerRadius^2)*rho1...
            + pi*(outerRadius^2-innerRadius^2)*rho2...
            +(4*(xmax-xmid)*(ymax-ymid)...
            - pi*outerRadius^2)*rho3;    % analytical normalization factor
        % Note: If rho is not already normalized over the domain of interest, the
        % normalization factor rhoNorm can be found by integrating rho over
        % the domain of interest.
        rho = rho / rhoNorm;    % normalized distribution
        
        % Generating Swarm Blob Function:
        G_i = @(r_i) exp( - ( (X-r_i(1)).^2 + (Y-r_i(2)).^2 ) / (2*delta^2) )...
            / (2*pi*delta^2);            	% gaussian blob for each robot
        rho_N = @(r) 0;
        for i = 1:N
            rho_N = @(r) rho_N(r) + G_i(r(2*i-1:2*i));
        end
        rho_N = @(r) rho_N(r)/N;            % normalized swarm blob function
        Drho = @(r) abs(rho_N(r) - rho);    % distribution difference
        e_N = @(r) sum(sum(Drho(r)))*dx*dy; % error metric
        
        % Initial Guess for Solver (Robots within Domain):
        x0 = rand(1,2*N);
        x0(1:2:end) = (xmax-4*delta)*x0(1:2:end) + (xmin + 2*delta);
        x0(2:2:end) = (ymax-4*delta)*x0(2:2:end) + (ymin + 2*delta);
        
        % Solver Parameters:
        f = e_N;            	% objective function
        A = [];                 % linear inequality constraint matrix
        b = [];                 % linear inequality constraint vector
        Aeq = [];               % linear equality constraint matrix
        beq = [];               % linear equality constraint vector
        lb = [xmin,ymin];       % lower bound on decision variables
        lb = repmat(lb,1,N);    % lower bound on all N robots
        ub = [xmax,ymax];       % upper bound on decision variables
        ub = repmat(ub,1,N);    % upper bound on all N robots
        c = [];                 % nonlinear constraints
        Fmax = 1000000;         % maximum functions calls by solver
        Imax = 1000000;         % maximum solver iterations
        
        % Optimization Solution:
        timeOfStart = datestr(now,30);  % start time of solver
        name = ['optimalRingIntermediateData',timeOfStart,'.txt'];
        file = fopen(name,'a');         % create txt file for intermediate data
        options = optimoptions('fmincon','MaxFunctionEvaluations',Fmax,...
            'Maxiterations',Imax,'OutputFcn',...
            @(x,optimValues,state)saveIteration(x,optimValues,state,file));
        [x,fval,exitflag,output] = fmincon(f,x0,A,b,Aeq,beq,...
            lb,ub,c,options);           % solver
        intermediateData = load(name);	% loads intermediate data to workspace
        
        % Figure Generation - Optimal Swarm Blob Function:
        figure();
        surf(X,Y,rho_N(x));
        view([1,1,1]);
        shading interp;
        xlim([xmin,xmax]);
        ylim([ymin,ymax]);
        ax = gca;
        ax.XTick = xmin:(xmax-xmin)/4:xmax;
        ax.YTick = xmin:(ymax-ymin)/7:ymax;
        pbaspect([(xmax-xmin) (ymax-ymin) (xmax-xmin)/2]);
        ttl = ['\begin{tabular}{c}Optimal Distribution for $N=$ ',...
            num2str(N),' Robots, $\delta = $ ',num2str(delta),' \\'...
            '$e_N^{\delta}=$ ',num2str(fval),'\end{tabular}'];
        title(ttl,'interpreter','latex');
        xlabel('$x$','interpreter','latex');
        ylabel('$y$','interpreter','latex');
        set(gca,'ticklabelinterpreter','latex');
        
        % Figure Generation - Solver Convergence Plot:
        figure();
        plot(intermediateData(:,1),intermediateData(:,2));
        ttl = ['\texttt{fmincon} Convergence for $N=$ ',...
            num2str(N),' Robots'];
        title(ttl,'interpreter','latex');
        xlabel('Iteration, $i$','interpreter','latex');
        ylabel('Lower Error Bound, $e_N^{\delta}$','interpreter','latex');
        set(gca,'ticklabelinterpreter','latex');
        
        % Saving Optimization Results:
        save(['optimalRingDistribution',timeOfStart,'.mat']);   % save workspace
    catch
        disp('Something went wrong. The catch block has been entered.');    % warning message
    end
    
    % Note: The following code stores the optimal error values.
    errors(k) = fval;       % append error to vector of stored values
    save('errors.mat','errors');
end