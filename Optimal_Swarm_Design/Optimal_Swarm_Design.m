clc;
clear;
close all;
rng shuffle;    % shuffle the random number generator seed for multistart on different computers

%%% Optimal Swarm Design Calculator

% This script computes the error metric minimum for a given number of robots N.
% The effective robot radius, delta, is included as a decision variable to
% find the optimal delta value corresponding to that N.

% To use this script, the user must enter the parameters noted by the USER
% INPUT sections below. These inputs include the swarm parameters, the
% computational parameters, and the target distribution.

% USER INPUT - Number of Robots:
N = [4,7,10];           % desired swarm sizes to optimize for

% Looping Over Different Swarm Sizes:
parfor j = 1:length(N)
    try
        disp(['Optimization for N = ',num2str(N(j)),' started.']);	% start message
        
        % Check for Continuation of Existing Optimization:
        name = ['optimizationDataN',num2str(N(j)),'.txt'];          % txt file for data storage
        if exist(name,'file') ~= 0
            prevData = load(name);
        end
        
        % USER INPUT - Swarm and Computational Parameters:
        xmin = 0;               % left boundary
        xmax = 4.8;             % right boundary
        ymin = 0;               % bottom boundary
        ymax = 7;               % top boundary
        nx = 200;               % number of grid points, x
        ny = 200;               % number of grid points, y
        delta0 = 0.1;           % initial guess for optimal delta
        deltaLb = 0;          	% lower bound on delta
        deltaUb = 0.7;       	% upper bound on delta
        
        % Computed Parameters:
        xmid = (xmin+xmax)/2;   % center, x
        ymid = (ymin+ymax)/2;   % center, y
        dx = (xmax-xmin)/nx;    % grid spacing, x
        dy = (ymax-ymin)/ny;    % grid spacing, y
        firstIter = 0;          % iteration starting number
        if N(j) > 200
            deltaUb = 0.15;     % tighten upper bound for large N
        elseif N(j) > 20 && N(j) < 200
            deltaUb = 0.45;     % tighten upper bound for medium N
        end
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
            - pi*outerRadius^2)*rho3;    % analytically derived normalization factor
        % Note: If rho is not already normalized over the domain of interest, the
        % normalization factor rhoNorm can be found by integrating rho over
        % the domain of interest.
        rho = rho / rhoNorm;    % normalized desired distribution
        
        % Generating Swarm Blob Function:
        G_i = @(r_i) exp( - ( (X-r_i(2)).^2 + (Y-r_i(3)).^2 ) / (2*r_i(1)^2) )...
            / (2*pi*r_i(1)^2);  % gaussian blob for each robot
        rho_N = @(r) 0;
        for i = 1:N(j)
            rho_N = @(r) rho_N(r) + G_i([r(1),r(2*i:2*i+1)]);
        end
        rho_N = @(r) rho_N(r)/N(j);       	% normalized swarm blob function
        Drho = @(r) abs(rho_N(r) - rho);    % distribution difference
        e_N = @(r) sum(sum(Drho(r)))*dx*dy; % error metric
        
        % Initial Guess for Solver (Robots within Domain):
        x0 = rand(1,2*N(j));
        x0(1:2:end) = xmax*x0(1:2:end);
        x0(2:2:end) = ymax*x0(2:2:end);
        x0 = [delta0,x0];
        
        % Initial Guess for Solver (Robots Defined by Previous Optimization):
        % Note: Use the following initial guess if continuing a previous
        % optimization.
        % x0 = prevData(end,3:end);       % set initial guess to last
        % firstIter = prevData(end,1);    % start iteration number from last
        % x0 = [delta0,x0];               % set initial delta guess to last
        
        % Solver Parameters:
        f = e_N;            	% objective function
        A = [];                 % linear inequality constraint matrix
        b = [];                 % linear inequality constraint vector
        Aeq = [];               % linear equality constraint matrix
        beq = [];               % linear equality constraint vector
        lb = [xmin,ymin];       % lower bound on decision variables
        lb = repmat(lb,1,N(j));	% lower bound on all N robots
        lb = [deltaLb,lb];      % lower bound for all decision variables
        ub = [xmax,ymax];       % upper bound on decision variables
        ub = repmat(ub,1,N(j));	% upper bound on all N robots
        ub = [deltaUb,ub];      % upper bound for all decision variables
        c = [];                 % nonlinear constraints
        Fmax = 1000000;         % maximum functions calls by solver
        Imax = 1000000;         % maximum solver iterations
        
        % Optimization Solution:
        tic;                            % start timing the solution
        file = fopen(name,'a');         % open txt file for intermediate data
        timeLim = 86400;                % time limit for fmincon call
        options = optimoptions('fmincon','MaxFunctionEvaluations',Fmax,...
            'Maxiterations',Imax,'OutputFcn',...
            @(x,optimValues,state)saveIteration(x,optimValues,state,file,timeLim,firstIter));
        [x,fval,exitflag,output] = fmincon(f,x0,A,b,Aeq,beq,...
            lb,ub,c,options);           % solver
        intermediateData = load(name);	% loads intermediate data to workspace
        timeElapsed = toc;              % store solution time
        disp(['Optimization for N = ',num2str(N(j)),...
            ' complete. Time elapsed: ',num2str(timeElapsed),' s.']);   % completion message
        
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
        pbaspect([(xmax-xmin),(ymax-ymin),(xmax-xmin)/2]);
        ttl = ['\begin{tabular}{c}Optimal Distribution for $N=$ ',...
            num2str(N(j)),' Robots, $\delta = $ ',num2str(x(1)),' \\'...
            '$e_N^{\delta}=$ ',num2str(fval),'\end{tabular}'];
        title(ttl,'interpreter','latex');
        xlabel('$x$','interpreter','latex');
        ylabel('$y$','interpreter','latex');
        set(gca,'ticklabelinterpreter','latex');

        % Figure Generation - Solver Convergence Plot:
        figure();           % intermediate solver data plot
        plot(intermediateData(:,1),intermediateData(:,2));
        ttl = ['\texttt{fmincon} Convergence for $N=$ ',...
            num2str(N(j)),' Robots, $\delta = $ ',num2str(x(1))];
        title(ttl,'interpreter','latex');
        xlabel('Iteration, $i$','interpreter','latex');
        ylabel('Error, $e_N^{\delta}$','interpreter','latex');
        set(gca,'ticklabelinterpreter','latex');
        
        % Saving Optimization Results:
        parsave(x(1),x,N(j),fval,X,Y,rho_N(x),...
            xmin,xmax,ymin,ymax,intermediateData,timeElapsed);     % save workspace
    catch
        disp('Something went wrong. The catch block has been entered.');    % warning message
    end
end

% Note: The following code stores the prescribed N values and their corresponding
% optimal errors into an array and saves it.
NvEdata = zeros(length(N),2);
for j = 1:length(N)
    intData = load(['optimizationDataN',num2str(N(j)),'.txt']); % load optimization data
    NvEdata(j,1) = N(j);                % store N
    NvEdata(j,2) = intData(end,2);      % store error
end
save('NvEdata.mat','NvEdata');          % save data
