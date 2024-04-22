% Main Script to run SEM 2D solver for the GS problem.
% 
% Assembles and solves the 2D SEM system for the problem
%           Lu(x,y) = f(x,y) in [a,b]x[c,d]
% with
%           Lu = -div(mu(x,y)*grad(u(x,y))) + sigma(x,y)*u
% and BC   u = g(x,y) on the boundary
% 
% The L operator is chosen to be the Grad-Shafranov operator
% divided by r:
%           mu(x,y) = 1/x,  sigma(x,y) = 0


%% Initialization
clearvars, clc, close all
addpath ./functions
addpath ./grad_shafranov
verbose = 0; % plots mesh and other stuff

input_file = 'input_data_JT60SA.mat';
in = load(input_file);

%% Problem definition
problem = 5;
switch problem
  case 1 % Homogeneous BC, for testing with FD poisson solver
    mu    = @(X,Y) ones(size(X));
    sigma = @(X,Y) zeros(size(X));
    f     = @(X,Y) X.^2 + Y.^2;
    g     = @(X,Y) zeros(size(X));
  case 5 % Grad-Shafranov problem with non-zero BC - TBD
    mu    = @(X,Y) min(1./X,1e10); 
%     mu    = @(X,Y) mu_fcn(X,Y); % mu_fcn is defined below
    sigma = @(X,Y) zeros(size(X));
    f     = @(X,Y) -RHS_GS(X,Y,in);
    g     = @(X,Y) BC_GS(X,Y,in);
end

%% Geometry definition
G.N    = 12;                        % order of the integration
G.MR   = 10;                        % no. of elements along y direction
G.MC   = 10;                        % no. of elements along x direction

G.xmin = 0;                        % domain definition (JT-60SA)
G.xmax = 10;            
G.ymin = -10; 
G.ymax = 10; 

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
[Xr,Yr,Ur] = refine_solution_SEM_2D(G,EE,un,gi,20,20,1); 

%% Benchmark with CNL
figure
subplot(121)
plot3(Xr,Yr,Ur,'ob')
hold on, axis equal
benchmark = load('CNL_solution.mat');
CNLsol = griddata(benchmark.Input_struct.p(1,:),benchmark.Input_struct.p(2,:),benchmark.x_np(1:size(benchmark.Input_struct.p,2)),Xr,Yr);
plot3(Xr,Yr,CNLsol,'.r');

subplot(122)
mesh(Xr,Yr,Ur-CNLsol)

%% Plot for GS
figure
subplot(121)

contourf(Xr,Yr,Ur,30);

subplot(122)
plot3(Xr,Yr,RHS_GS(Xr,Yr,in),'.r')


%% Check against finite-differences Poisson solver
if problem==1
  addpath ../poisson_test
  figure
  hold on
  fh = testFdPoisson(G.xmin,G.xmax,f,g);
  rmpath ../poisson_test
  hold on
  figure(fh)
  plot3(xn,yn,un,'.r','linewidth',2,'markersize',20)
  plot3(xb,yb,(un(gb)),'ob','linewidth',2)
end





function out = mu_fcn(X,Y)
    out = 1./X;
    out(isnan(out)) = 0;
    out(isinf(out)) = 0;
end
