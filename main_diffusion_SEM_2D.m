error('solutions are different if parameters change! CHECK!')



%% Initialization
clearvars, clc, close all
addpath ./functions
verbose = 0; % plots mesh and other stuff

%% Problem definition
mu    = @(X,Y) ones(size(X));
sigma = @(X,Y) zeros(size(X));
f     = @(X,Y) ones(size(X));
g     = @(X,Y) zeros(size(X));
U0fun = @(X,Y) gauss2D(X,Y,0.3,0.5);

%% Time interval
TT = 3;
tspan = [0 TT];

%% Geometry definition
G.N    = 2;                         % order of the integration
G.MR   = 4;                         % no. of elements along y direction
G.MC   = 4;                         % no. of elements along x direction

G.xmin = -1;                        % domain definition
G.xmax = 1;
G.ymin = G.xmin;
G.ymax = G.xmax;

G.NN = (G.MR*G.N+1)*(G.MC*G.N+1);   % total no. of nodes
G.NB = (G.MR*G.N)*2+(G.MC*G.N)*2-4; % no. of boundary nodes

%% Solve problem
[EE,xn,yn,gi,xb,yb,gb,xm,ym,gm] = initmeshSEM(G,verbose);
gin = setdiff(gi,gb); % inner points

AA = SEM_2D_operator(G,EE,g,mu,sigma);
U0 = U0fun(xn,yn); % initial condition

% ops = odeset('MaxStep',1e-3,'RelTol',1e-2);
% [t,Un] = ode23tb(@(t,U) odefun_diffusion2D_SEM(t,U,AA,g,gin,gb,xb,yb), tspan, U0, ops);

ops = odeset('MaxStep',1e-1,'RelTol',1e-2);
[t,Un] = ode45(@(t,U) odefun_diffusion2D_SEM(t,U,AA,g,gin,gb,xb,yb), tspan, U0, ops);

%% Plot
figure()
plot3(xn,yn,U0,'.r','Markersize',20);
aa = [min(xn), max(xn), min(yn), max(yn), 0, max(max(U0))];
axis(aa)

subplot(221)
[x,y,~] = reshape_solution_SEM_2D(xn,yn,U0,G);
surf(x,y,zeros(size(x)))
axis(aa(1:4));

nstep = 10;
for i = 1 : nstep : length(t)
  
  [x,y,u] = reshape_solution_SEM_2D(xn,yn,Un(i,:),G);
  
  subplot(222)
  contour(x,y,u)
  
  subplot(223)
  plot3(xn,yn,Un(i,:),'.r','Markersize',10);
  hold on
  mesh(x,y,u);
  axis(aa);
  hold off
  
  subplot(224)
  plot3(xn,yn,Un(i,:),'.r','Markersize',10);
  hold on
  refine_solution_SEM_2D(G,EE,Un(i,:),gi,20,20,1);
  title(['t = ' num2str(t(i))])
  axis(aa);
  hold off
  
  pause(1e-2);
  
end

%==========================================================================

%% Main ODE function
function dUdt = odefun_diffusion2D_SEM(~,U,AA,g,gin,gb,xb,yb)
dUdt(gin,1) = -AA(gin,gin)*U(gin);
dUdt(gb,1)  = g(xb,yb);
end