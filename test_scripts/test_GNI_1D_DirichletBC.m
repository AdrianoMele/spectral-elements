%% Test G-NI 1D
% Solve      Lu = -(mu')' + su = f in [-1,1]
% with BC    u(-1) = u(1) = 0
% u >= u0 > 0, s >= 0 to have a bilinear form which is coercive
% (Sample problem in Quarteroni et al. p.236)

clearvars
close all
clc

addpath functions

%% Problem definition

problem = 1;

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
    f     = @(x) x.^2 - 3;
    sol   = @(x) x.^2 - 1;
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

% Domain
dx = 1e-4;
x  = -1:dx:1;

% Order of the integration
N = 7; 

% GLL quadrature nodes
xk = find_GLL_nodes(N);

% Characteristic polynomials for GLL nodes
Lk = legendre_poly(N,xk);
Psi = GLL_poly(N,x,xk,Lk);

% Integration coefficients (analytical)
ak = 2/(N*(N+1)) ./ (Lk'.^2); 
% ak = sum(Psi,2)*dx;           % Numerical


%% Compute RHS and stiffness matrix
Ppk = GLL_poly_derivatives(N,xk,xk,Lk);

% RHS
F = ak.*f(xk)';

% Stiffness matrix (u=0 at boundary)
A = zeros(N+1);
for j = 1:N+1
  for i = 1:N+1
    A(j,i) = scalprod(mu(xk).*Ppk(j,:), Ppk(i,:), [], ak);
  end  
  A(j,j) = A(j,j) + sigma(xk(j))*ak(j);
end
% A = A + diag(sigma(xk).*ak);

% BC
PpN = Ppk(:,end); 
Pp0 = Ppk(:,1);
for j = [1, N+1]
  for i = 1:N+1
    A(j,i) = A(j,i) - mu(1)*PpN(i)*Psi(j,end) + mu(-1)*Pp0(i)*Psi(j,1);
  end
end
F = F - A(:,1)*u1 - A(:,N+1)*u2; 
A = A(:,2:N);

% Find solution
un = A\F;
un = [u1; un; u2];

% Plot interpolant
u = x*0;
for i = 1:N+1
  u = un(i)*Psi(i,:) + u;
end

plot(x,u,'r', xk,un,'^k', x,sol(x),'--b'), grid on
title(sprintf('RMSE: %d', rms(u-sol(x))))

rmpath functions


