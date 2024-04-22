%% Test G-NI 1D with mixed boundary conditions
% Solve      Lu = -(mu')' + su = f in [-1,1]
% with BC    u(-1) = 0, u'(1) = 1
% u >= u0 > 0, s >= 0 to have a bilinear form which is coercive
% (Sample problem in Quarteroni et al. p.236)

clearvars
close all
clc

%% Problem definition

problem = 4;

% Problem parameters (depending on x), RHS, analytic solution
switch problem
  case 4 % Quarteroni, sample problem p.254 (diffusion-reaction)
    mu    = @(x) 1+x.^2;
    sigma = @(x) cos(x.^2);
    f     = @(x) -(x.^2);
    sol   = @(x) NaN;
    u1 = 0; up2 = 1; % mixed BC    
end


%% Initialization

% Domain
dx = 1e-4;
x  = -1:dx:1;

% Order of the integration
N = 9; 

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
F = ak.*f(xk)' + mu(1)*Psi(:,end);

% Stiffness matrix (u=0 at boundary)
A = zeros(N+1);
for i = 1:N+1
  for j = 1:N+1
    A(i,j) = scalprod(mu(xk).*Ppk(j,:), Ppk(i,:), [], ak);
  end  
  A(i,i) = A(i,i) + sigma(xk(i))*ak(i);
end
% A = A + diag(sigma(xk).*ak);

% for i = 2:N+1
%   A(i,i) = A(i,i) + sigma(xk(i))*ak(i);
% end


% Find solution
un = A(2:N+1,2:N+1)\F(2:N+1);

% Find interpolant
u = x*0;
for i = 1:N
  u = un(i)*Psi(i+1,:) + u;
end


%% Show results

figure
plot(x,u,'r', xk,[0; un],'^k', x,sol(x),'--b'), grid on
title(sprintf('RMSE: %d', rms(u-sol(x))))

% frec = - diff(mu(x(1:end-1).*diff(u)./diff(x)))./diff(x(1:end-1)) + sigma(x(2:end-1)).*u(2:end-1);
% figure
% plot(x,f(x),x(2:end-1),frec,'--r')


