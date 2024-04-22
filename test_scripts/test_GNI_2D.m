%% Test G-NI 2D
% Solve      Lu = -div(m*grad(u)) + s*u = f in [-1,1]x[-1,1]
% with BC    u = 0 on the boundary
% mu >= mu0 > 0, s >= 0 to have a bilinear form which is coercive
% (Sample problem in Quarteroni et al. p.245)

clearvars
close all
clc

addpath ../functions

%% Problem definition
problem = 1;
switch problem
  case 1
    mu    = @(X,Y) ones(size(X));
    sigma = @(X,Y) zeros(size(X));
    f     = @(X,Y) X.^2 + Y.^2;
    g     = @(X,Y) X*0;
    sol   = @(X,Y) NaN;  % Check consistency with finite differences Poisson solver
end

%% Initialization

% Reference domain
dx = 5e-2;    dy = dx;
x  = -1:dx:1; y = x;

% Order of the integration
N  = 5;

[xk,ak,Lk,XK,YK,ajk] = init_ref_element(N);


%% Compute RHS and stiffness matrix
Psi2D           = GLL_poly_2D(N,x,y,xk,Lk);
[Ppkx2D,Ppky2D] = GLL_poly_2D_gradient(N,xk,xk,xk,Lk);

F  = zeros(N+1,N+1);
A1 = zeros(N+1,N+1,N+1,N+1);
A2 = zeros(N+1,N+1,N+1,N+1);

% Stiffness matrix and RHS
% The product takes the form
% A[l,r x i,j] u[i,j] = F[l,r]
%
%#ok<*ALIGN>
for l = 1:N+1, for r = 1:N+1 
  for i = 1:N+1, for j = 1:N+1  
    
    for m = 1:N+1, for n = 1:N+1 
      A1(l,r,i,j) = A1(l,r,i,j) + mu(xk(m),xk(n)).*ajk(m,n).*(Ppkx2D(i,j,m,n).*Ppkx2D(l,r,m,n) + Ppky2D(i,j,m,n).*Ppky2D(l,r,m,n));
      % A2(l,r,i,j) = A2(l,r,i,j) + sigma(xk(m),xk(n)) * kronecker(i,l) * kronecker(j,r);           
    end, end
    
    A2(l,r,l,r) = A2(l,r,l,r) + sigma(xk(l),xk(r));
    end, end

  F(l,r) = f(xk(l),xk(r))*ajk(l,r);
end, end

A = A1+A2;


%% Assemble problem
FF = zeros((N-1)^2,1); % Remove boundary
AA = zeros((N-1)^2);
for l = 1:N-1
  for r = 1:N-1
    FF((l-1)*(N-1)+r) = F(l+1,r+1);   
    for i = 1:N-1
      for j = 1:N-1
        AA((l-1)*(N-1)+r,(i-1)*(N-1)+j) = A(l+1,r+1,i+1,j+1);
      end
    end   
  end
end

% Solve
uij = AA\FF;


%% Plot comparison with finite differences solution
% [XK,YK] = meshgrid(xk(2:N),xk(2:N));
Uij = zeros(N-1);
for i = 1:N-1
  for j = 1:N-1
    Uij(i,j) = uij((i-1)*(N-1)+j);
  end
end

% figure
% surf(XK,YK,Uij)

[X,Y] = meshgrid(x,y);
U = zeros(length(x));
for i = 1:N-1
  for j = 1:N-1
    U = U + Uij(i,j).*squeeze(Psi2D(i+1,j+1,:,:));
  end
end

addpath poisson_test
fh = testFdPoisson(-1,1,f,g);
rmpath poisson_test
hold on
plot3(X,Y,U,'.r');

rmpath functions
