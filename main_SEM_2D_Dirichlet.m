%% Initialization
clearvars, clc, close all
addpath ./functions
verbose = 0; % plots mesh and other stuff

%% Problem definition
problem = 4;
switch problem
  case 1 % Homogeneous BC
    mu    = @(X,Y) ones(size(X));
    sigma = @(X,Y) zeros(size(X));
    f     = @(X,Y) X.^2 + Y.^2;
    g     = @(X,Y) zeros(size(X));
  case 2 % non-zero BC
    mu    = @(X,Y) ones(size(X));
    sigma = @(X,Y) zeros(size(X));
    f     = @(X,Y) -5*pi^2*sin(pi*X).*cos(2*pi*Y);
    g     = @(X,Y) -sin(pi*X).*cos(2*pi*Y);
  case 3 % Same as 2, but swap coordinates to check BC
    mu    = @(X,Y) ones(size(X));
    sigma = @(X,Y) zeros(size(X));
    f     = @(X,Y) -5*pi^2*sin(pi*Y).*cos(2*pi*X);
    g     = @(X,Y) -sin(pi*Y).*cos(2*pi*X);
  case 4 % Homogeneous BC
    mu    = @(X,Y) ones(size(X));
    sigma = @(X,Y) zeros(size(X));
    f     = @(X,Y) ones(size(X));
    g     = @(X,Y) zeros(size(X));
end

%% Geometry definition
G.N    = 5;                         % order of the integration
G.MR   = 2;                         % no. of elements along y direction
G.MC   = 2;                         % no. of elements along x direction

G.xmin = -1;                        % domain definition
G.xmax = 1;            
G.ymin = G.xmin; 
G.ymax = G.xmax; 

G.NN = (G.MR*G.N+1)*(G.MC*G.N+1);   % total no. of nodes 
G.NB = (G.MR*G.N)*2+(G.MC*G.N)*2-4; % no. of boundary nodes


%% Solve problem
[EE,xn,yn,gi,xb,yb,gb,xm,ym,gm] = initmeshSEM(G,verbose);
un = solve_SEM_2D(G,EE,g,f,mu,sigma,xb,yb,gb);
[x,y,u] = reshape_solution_SEM_2D(xn,yn,un,G);


%% Plot solution
figure();
subplot(121)
hold on
plot3(xn,yn,un,'.r','linewidth',2,'markersize',20)
plot3(xb,yb,(un(gb)),'ob','linewidth',2)
mesh(x,y,u);
subplot(122)
refine_solution_SEM_2D(G,EE,un,gi,40,40,1); % 100 points per element



%% Check against finite-differences Poisson solver
addpath ../poisson_test
figure
hold on
fh = testFdPoisson(G.xmin,G.xmax,f,g);
rmpath ../poisson_test
hold on
figure(fh)
plot3(xn,yn,un,'.r','linewidth',2,'markersize',20)
plot3(xb,yb,(un(gb)),'ob','linewidth',2)


