%% Test G-NI 2D
% Solve      Lu = -div(m*grad(u)) + s*u = f in [-1,1]x[-1,1]
% with BC    u = g on the boundary
% mu >= mu0 > 0, s >= 0 to have a bilinear form which is coercive
% (Sample problem in Quarteroni et al. p.245)

clearvars
close all
clc
addpath functions

verbose = 0;

%% Problem definition
problem = 1;
switch problem
  case 1
    mu    = @(X,Y) ones(size(X));
    sigma = @(X,Y) zeros(size(X));
    f     = @(X,Y) X.^2 + Y.^2;
    g     = @(X,Y) zeros(size(X));
  case 2
    mu    = @(X,Y) ones(size(X));
    sigma = @(X,Y) zeros(size(X));
    f     = @(X,Y) -5*pi^2*sin(pi*X).*cos(2*pi*Y);
    g     = @(X,Y) -sin(pi*X).*cos(2*pi*Y);
%   case 3 % not working
%     mu    = @(X,Y) ones(size(X));
%     sigma = @(X,Y) zeros(size(X));
%     f     = @(X,Y) -2*pi^2*sin(pi*X^2).*cos(pi*Y);
%     g     = @(X,Y) -sin(2*pi*X).*cos(2*pi*Y);
end


%% Initialization

% Reference domain
dx = 5e-2;    dy = dx;
x  = -1:dx:1; y = x;

% Order of the integration
N = 2;

[xk,ak,Lk,XK,YK,ajk] = init_ref_element(N);

%% Actual domain
a = 0; b = 2;
xve = [a b b a];
yve = [a a b b];

[XEK,YEK,J,M] = affine_transform(xve,yve,XK,YK);

ajk = ajk*M/4;


%% Compute RHS and stiffness matrix
Psi2D           = GLL_poly_2D(N,x,y,xk,Lk);
[Ppkx2D,Ppky2D] = GLL_poly_2D_gradient(N,xk,xk,xk,Lk);

F  = zeros(N+1,N+1);
A1 = zeros(N+1,N+1,N+1,N+1);
A2 = zeros(N+1,N+1,N+1,N+1);
A3 = zeros(N+1,N+1,N+1,N+1);

% Stiffness matrix and RHS
% The product takes the form
% A[l,r x i,j] u[i,j] = F[l,r]
%
%#ok<*ALIGN>
for l = 1:N+1, for r = 1:N+1    
  for i = 1:N+1, for j = 1:N+1     
    for m = 1:N+1, for n = 1:N+1 
      A1(l,r,i,j) = A1(l,r,i,j) + mu(XEK(n,m),YEK(n,m)).*ajk(n,m).*(Ppkx2D(i,j,n,m).*Ppkx2D(l,r,n,m) + Ppky2D(i,j,n,m).*Ppky2D(l,r,n,m))/J(n,m);
      % A2(l,r,i,j) = A2(l,r,i,j) + sigma(xk(m),xk(n)) * kronecker(i,l) * kronecker(j,r);           
    end, end   
    A2(l,r,l,r) = A2(l,r,l,r) + sigma(XEK(l,r),YEK(l,r));
    end, end
  F(l,r) = f(XEK(l,r),YEK(l,r))*ajk(l,r);
end, end

% BC - check weights!
for i = 1:N+1, for j = 1:N+1
    for m = 1:N+1 % This m should be a l
        
    % Bottom
    A3(1,m,i,j)   = A3(1,m,i,j)   + mu(XEK(1,m),YEK(1,m)).*ak(m).*Ppky2D(i,j,1,m);   
    % Top
    A3(N+1,m,i,j) = A3(N+1,m,i,j) - mu(XEK(N+1,m),YEK(N+1,m)).*ak(m).*Ppky2D(i,j,N+1,m);
    % Left
    A3(m,1,i,j)   = A3(m,1,i,j)   + mu(XEK(m,1),YEK(m,l)).*ak(m).*Ppkx2D(i,j,m,1);
    % Right
    A3(m,N+1,i,j) = A3(m,N+1,i,j) - mu(XEK(m,N+1),YEK(m,N+1)).*ak(m).*Ppkx2D(i,j,m,N+1);
    
   end, end  
end

A = A1+A2+A3;

xek = zeros((N+1)^2,1); yek = xek;
% idx = [];
for i = 1:N+1, for j = 1:N+1
%     idx = [idx; [(i-1)*(N+1)+j i j] ];
    xek((i-1)*(N+1)+j) = XEK(i,j);
    yek((i-1)*(N+1)+j) = YEK(i,j);
end, end

% Find indexes of the boundary nodes
idxl = ((1:N+1)-1)*(N+1)+1;
idxr = ((1:N+1)-1)*(N+1)+(N+1);
idxb = (1-1)*(N+1) + (1:N+1);
idxt = (N+1-1)*(N+1) + (1:N+1);
idxbound = unique([idxt idxr idxb idxl]);
idxin    = setdiff(1:(N+1)^2,idxbound); 

xb = xek(idxbound); yb = yek(idxbound);
ub = g(xb,yb);
% plot(xb,yb,'ob');


%% Assemble problem
FF = zeros((N+1)^2,1); 
AA = zeros((N+1)^2,(N-1)^2);

for l = 1:N+1, for r = 1:N+1
  FF((l-1)*(N+1)+r) = F(l,r);
  
  for i = 1:N+1, for j = 1:N+1
    AA((l-1)*(N+1)+r,(i-1)*(N+1)+j) = A(l,r,i,j);
  
  end, end
end, end

FF = FF - AA(:,idxbound)*ub;
AA = AA(:,idxin);


%% Solve
uij = zeros((N+1)^2,1);
uij(idxin) = AA\FF;
uij(idxbound) = ub;


%% Plot comparison with finite differences solution
% [XK,YK] = meshgrid(xk(2:N),xk(2:N));
Uij = zeros(N-1);
for i = 1:N+1
  for j = 1:N+1
    Uij(i,j) = uij((i-1)*(N+1)+j);
  end
end
% Uij = [ubb; ubl(2:N) Uij ubr(2:N); ubt];

[X,Y] = meshgrid(x,y);
[XE,YE] = affine_transform(xve,yve,X,Y);
U = zeros(length(x));
for i = 1:N+1
  for j = 1:N+1
    U = U + Uij(i,j).*squeeze(Psi2D(i,j,:,:));
  end
end

addpath poisson_test
fh = testFdPoisson(a,b,f,g);
rmpath poisson_test
hold on
plot3(XE,YE,U,'.r','linewidth',2);
plot3(XEK,YEK,Uij,'ob')
plot3(XEK,YEK,g(XEK,YEK),'+k') %BC

if verbose
  for i = 1 : N+1
    for j = 1 : N+1
      clf
      mesh(XE,YE,squeeze(Psi2D(i,j,:,:)))
      hold on
      plot3(XEK(i,j),YEK(i,j),1,'ob')
      plot3(xek(j),xek(i),1,'+r')
      pause
    end
  end
end

rmpath functions