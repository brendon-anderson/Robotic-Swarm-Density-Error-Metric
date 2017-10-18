clear;
close all;
clc;

%%% Error Metric Probability Density Function Calculator

% This script computes the CDF of the error for robots over a target
% distribution using monte carlo integration. The code then curve fits the
% CDF, then differentiates the curve fit expression to obtain a PDF of the
% error.

% To use this script, the user must enter the parameters noted by the USER
% INPUT sections below. These inputs include the swarm parameters, the
% computational parameters, and the target distribution in the variable
% definitions below.

% USER INPUT - Swarm and Computational Parameters:
nRobots = 200;          % number of robots in swarm
delta = 0.2;            % effective robot radius
nMC = 500;              % number of  monte carlo evaluation points
nx = 200;               % number of spatial grid points, x
ny = 200;               % number of spatial grid points, y
nE = 50;                % number of error metric evaluation points
nt = 1000;              % number of grid points for CDF/PDF plots
x0 = [0.5,20,0.77,0.5];	% initial coefficient guess for CDF curve fit of form: CDF = x(1)*erf(x(2)*(error-x(3)))+x(4)
xmin = 0;               % left spatial boundary
xmax = 4.8;             % right spatial boundary
ymin = 0;               % bottom spatial boundary
ymax = 7;               % top spatial boundary

% Computed Parameters:
nDim = 2*nRobots;           % number of vector dimensions
xmid = (xmin+xmax)/2;       % center, x
ymid = (ymin+ymax)/2;       % center, y
dx = (xmax-xmin)/nx;        % grid spacing, x
dy = (ymax-ymin)/ny;        % grid spacing, y
X = xmin:dx:xmax;           % x grid
Y = ymin:dy:ymax;           % y grid
[X,Y] = meshgrid(X,Y);      % grid
Xtens = repmat(X,1,1,nMC);  % tensor of X grids
Ytens = repmat(Y,1,1,nMC);  % tensor of Y grids

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
    - pi*outerRadius^2)*rho3;       % analytically derived normalization factor
% Note: If rho is not already normalized over the domain of interest, the
% normalization factor rhoNorm can be found by integrating rho over
% the domain of interest.
rho = rho / rhoNorm;                % normalizing target distribution
rho_i = @(r_i) (rho2-rho1)*double(((r_i(1)-xmid).^2 ...
    + (r_i(2)-ymid).^2 < outerRadius^2) & ...
    ((r_i(1)-xmid).^2 + (r_i(2)-ymid).^2 > innerRadius^2))...
    + rho1;                         % functional form of target distribution
% Note: rho_i takes a single robot position vector, r_i = [x_i,y_i], and
% computes the local target distribution value at that point.

% Generating Swarm Blob Function:
G_i = @(r_i) exp( - ( (bsxfun(@minus,Xtens,r_i(:,1,:))).^2 ...
    + (bsxfun(@minus,Ytens,r_i(:,2,:))).^2 ) / (2*delta^2) )...
    / (2*pi*delta^2);	% gaussian blob for each robot
rho_N = @(r) 0;
for i = 1:nRobots
    rho_N = @(r) rho_N(r) + G_i(r(:,2*i-1:2*i,:));
end
rho_N = @(r) rho_N(r)/nRobots;      % normalized swarm blob function
rho = repmat(rho,1,1,nMC);          % tensor of desired distributions
Drho = @(r) abs(rho_N(r) - rho);    % distribution difference
e_N = @(r) sum(sum(Drho(r)))*dx*dy; % error metric

% Placing Robots According to Target Density:
X0 = zeros(1,nDim,nMC);            	% robot positions vector
rhoX0 = zeros(1,nRobots,nMC);     	% scalar field values for each robot position
for j = 1:nMC
    a = zeros(nRobots,1);           % x vector
    b = zeros(nRobots,1);           % y vector
    for i = 1:nRobots
        x1=rand(1)*xmax;            % random x value
        y1=rand(1)*ymax;            % random y value
        z1=rand(1)*36;              % random value in [0,36]
        while z1 > rho_i([x1,y1])   % check if value exceeds local scalar field
            x1=rand(1)*xmax;        % replace x until it doesn't
            y1=rand(1)*ymax;        % replace y until it doesn't
            z1=rand(1)*36;          % update value
        end
        a(i)=x1;                    % store x coordinate
        b(i)=y1;                    % store y coordinate
        rhoX0(1,i,j) = rho_i([x1,y1])/rhoNorm;    % store scalar field values
    end
    X0(1,1:2:end,j) = a;            % x values
    X0(1,2:2:end,j) = b;            % y values
end

% Monte Carlo Computations:
e = e_N(X0);                        % compute errors for robot configurations
eExpected = sum(e)/nMC;             % compute mean (expected) error value
% NOTE: ee represents the error values over which the PDF is computed
eemid = eExpected - mod(eExpected*100,5)/100;   % middle of error metric grid
eemin = eemid - eemid/4;         	% min error metric value to compute for
eemax = eemid + eemid/3;            % max error metric value to compute for
dE = (eemax-eemin)/nE;              % error metric grid space
eeX = linspace(eemin,eemax,nE);     % error metric grid
eeX = repmat(eeX,1,1,nMC);          % error metric grid for N configurations
e = repmat(e,1,nE,1);
fX = zeros(1,nE,nMC);               % PDF values
fX(e <= eeX) = 1;                   % configurations in set
fX = permute(fX,[3,2,1]);           % move dimension three into rows
I = sum(fX)/nMC;                    % monte carlo mean
eeX = eeX(:,:,1);                	% restore ee domain

% CDF Curve Fitting:
fun = @(x,eeX) x(1)*erf(x(2)*(eeX-x(3)))+x(4);  % general erf form for CDF
x = lsqcurvefit(fun,x0,eeX,I);                  % solve for coefficients
t = linspace(eeX(1),eeX(end),nt);            	% grid for CDF/PDF plots
CDFfit = fun(x,t);                              % CDF curve fit
PDFfit = 2*x(1)*x(2)*exp(-(x(2)^2)*(t-x(3)).^2)/sqrt(pi);   % derivative of CDF curve fit

% Generated CDF and PDF Plots:
figure();
plot(eeX,I,t,CDFfit);       % plot computed and curve fit CDF
xlim([eemin,eemax]);
ttl = ['CDF for $N = $ ',num2str(nRobots),' Robots, $\delta = $ ',num2str(delta)];
title(ttl,'interpreter','latex');
xlabel('Error Metric, $e_N^{\delta}$','interpreter','latex');
ylabel('Cumulative Probability, $\mathcal{P}_e$','interpreter','latex');
legend({'Numerical Approximation','Curve Fit'},'interpreter','latex',...
    'location','northwest');
set(gca,'ticklabelinterpreter','latex');

figure();
plot(t,PDFfit);             % plot curve fit PDF
xlim([eemin,eemax]);
ttl = ['PDF for $N = $ ',num2str(nRobots),' Robots, $\delta = $ ',num2str(delta)];
title(ttl,'interpreter','latex');
xlabel('Error Metric, $e_N^{\delta}$','interpreter','latex');
ylabel('Probability Density, $\rho_e$','interpreter','latex');
set(gca,'ticklabelinterpreter','latex');