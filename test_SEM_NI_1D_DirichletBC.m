%% Test SEM-NI 1D
% Solve      Lu = -(mu')' + su = f in [-1,1]
% with BC    u(-1) = u(1) = 0
% u >= u0 > 0, s >= 0 to have a bilinear form which is coercive
% (Sample problem in Quarteroni et al. p.236)

clearvars
close all
% clc

addpath functions

verbose = 0;

%% Problem definition

problem = 2;

% Problem parameters (depending on x), RHS, analytic solution
switch problem
  case 1
    mu    = @(x) ones(size(x));
    sigma = @(x) zeros(size(x));
    f     = @(x) x.^2;
    sol   = @(x) -x.^4/12 + 1/12;
    u1 = 0; u2 = 0;
    
  case 2    
    m     = 2;
    mu    = @(x) m*ones(size(x));
    sigma = @(x) zeros(size(x));
    f     = @(x) x;
    sol   = @(x) -x.^3/(6*m) + x/(6*m);
    u1 = 0; u2 = 0;
  
  case 3
    m     = 1;
    s     = 1;
    mu    = @(x) m*ones(size(x));
    sigma = @(x) s*ones(size(x));
    f     = @(x) -(x.^2 - 3);
    sol   = @(x) -(x.^2 - 1);
    u1 = 0; u2 = 0;
    
  case 4
    m = 1;
    s = 2;
    mu    = @(x) m*ones(size(x));
    sigma = @(x) s*ones(size(x));
    f     = @(x) s*(x+1);
    sol   = @(x) x+1;
    u1 = 0; u2 = 2;
    
  case 5
    m = 1;
    s = 2;
    mu    = @(x) m*ones(size(x));
    sigma = @(x) s*ones(size(x));
    f     = @(x) -6*x + s*(x.^3+1);
    sol   = @(x) x.^3+1;
    u1 = 0; u2 = 2;
end


%% Initialization

%================== Reference element
np = 2000;
x = linspace(-1,1,np);
dx = x(2)-x(1);

% Order of the integration
N = 1;

[xk,ak,Lk] = init_ref_element(N);
Psi = GLL_poly(N,x,xk,Lk);

% % GLL quadrature nodes
% xk = find_GLL_nodes(N);
% 
% % Characteristic polynomials for GLL nodes
% Lk = legendre_poly(N,xk);
% Psi = GLL_poly(N,x,xk,Lk);
% 
% % Integration coefficients
% ak = 2/(N*(N+1)) ./ (Lk.^2); % Analytical
% % ak = sum(Psi,2)*dx;  % Numerical




%================= Move to finite elements

% Number of elements
M = 30;
xm = linspace(x(1),x(end),M+1);
% xm = [-1 2*sort(rand(1,M-1)-0.5) 1]; % random uneven spacing

% Nodes and weights for each element
tk  = zeros(M,N+1);
dk  = zeros(M,N+1);
xx  = zeros(M,length(x));
ddx = zeros(M,1);
for m = 1 : M
  a = xm(m);
  b = xm(m+1);
  
  tk(m,:) = ((b-a)*xk + (a+b))/2;    % new nodes
  dk(m,:) = (b-a)*ak/2;              % new coefficients
  
  xx(m,:)  = ((b-a)*x  + (a+b))/2;   % rescale x coordinate
  ddx(m,1) = x(2)-x(1);              % adapt x-step
end

% Plot subintervals and nodes
if verbose
  figure
  plot(xm,0*xm,'or',  tk,0*tk,'+k'), grid on, hold on
  plot(xx(1,:),Psi,'b')
end


%% Base functions

NP = 1 + (np-1)*M ;     % Total number of nodes
X = zeros(1,NP);
P = zeros(M*N+1,NP);
for m = 1 : M
  
  idxr = (m-1)*N + 1 : m*N + 1;
  idxc = (m-1)*(np-1) + 1 : m*(np-1)+1;
  
  X(idxc) = xx(m,:);    % Nodes associated to each element
  P(idxr,idxc) = Psi;   % Base functions extended to all nodes
end


% Plot base functions
if verbose
  figure
  for i = 1 : size(P,1)
    plot(X,P(i,:),tk,0*tk,'+k',xm,0*xm,'or'), grid on
    pause
  end
end

%% Compute RHS and stiffness matrix

% Derivatives of GLL polynomials (used in the stiffness matrix)
Ppk  = GLL_poly_derivatives(N,xk,xk,Lk);
dxm2 = diff(xm)/2;

% Compute matrices for each subinterval, then assemble
A = zeros(M*N+1);   % Stiffness matrix
F = zeros(M*N+1,1); % RHS
T = zeros(M*N+1,1); % nodes
for m = 1 : M
  idx = (m-1)*N+1 : m*N+1; % index of the nodes in the m-th interval
  
  % Nodes
  T(idx) = tk(m,:); 
  
  % Stiffness matrix
  Am = zeros(N+1);
  for j = 1:N+1
    for i = 1:N+1
      Am(j,i) = Am(j,i) + scalprod(mu(tk(m,:)).*Ppk(i,:)/dxm2(m), Ppk(j,:)/dxm2(m), [], dk(m,:));
    end
    Am(j,j) = Am(j,j) + sigma(tk(m,j))*dk(m,j);
  end
  A(idx,idx) = A(idx,idx) + Am;
  
  % Rhs
  fkm = f(tk(m,:)).*dk(m,:);
 F(idx) = F(idx) + fkm';
  
%   spy(A), pause
end

% BC
PpN = Ppk(:,end)/dxm2(M);
Pp0 = Ppk(:,1)  /dxm2(1);

A(1,1:N+1)             = A(1,1:N+1)             + mu(xm(1))*Pp0';
A(end,(M-1)*N+1:M*N+1) = A(end,(M-1)*N+1:M*N+1) - mu(xm(M))*PpN';

F = F - A(:,1)*u1 - A(:,end)*u2; 
A = A(:,2:end-1);

%% Compute solution
un = A\F;
un = [u1;un;u2];

% Compute interpolant
u = X*0;
for i = 1:size(P,1)
  u = un(i)*P(i,:) + u;
end


% Plot results
plot(T,un,'o',    X,u,'r',    X,sol(X),'--b',    xm,0*xm,'^b')
grid on
title(sprintf('RMSE: %d', rms(u-sol(X))))

rmpath functions


