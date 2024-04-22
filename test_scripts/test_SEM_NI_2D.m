%% Test SEM-NI 2D
% Solve      Lu = -div(m(x,y)*grad(u(x,y))) + s(x,y)*u(x,y) = f in [-1,1]x[-1,1]
% with BC    u = g(x,y) on the boundary
% u >= u0 > 0, s >= 0 to have a bilinear form which is coercive
% (Sample problem in Quarteroni et al. p.245)

clearvars
close all
clc
addpath ../functions
verbose = 1;

%% Problem definition
problem = 3;
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

%% Initialization

% Reference element =======================================================
dx = 5e-2;    dy = dx;
x  = -1:dx:1; y  = x;

% Order of the integration (after ~10 it becomes very slow)
G.N = 10;
[xk,ak,Lk,XK,YK,AK] = init_ref_element(G.N);

% GLL 2D polynomials over reference element (Nxk x Nyx x Nx x Ny)
Psi2D = GLL_poly_2D(G.N,xk,xk,xk,Lk);
[Ppkx2D,Ppky2D,XG,YG] = GLL_poly_2D_gradient(G.N,xk,xk,xk,Lk);

% Move to finite elements =================================================
G.xmin = -1; G.xmax = 1; G.ymin = G.xmin; G.ymax = G.xmax; % domain definition
G.MR = 2; % no. of columns of the grid
G.MC = 1; % no. of rows

% Total no. of nodes and boundary nodes
G.NN = (G.MR*G.N+1)*(G.MC*G.N+1);
G.NB = (G.MR*G.N)*2 + (G.MC*G.N)*2 - 4;

% Generate mesh
[EE,xn,yn,gi,xb,yb,gb,xm,ym,gm] = initmeshSEM(G,verbose);

%% Build stiffness matrix and RHS

% for test
FFF = zeros(G.NN,1);
AAA = zeros(G.NN);




% Stiffness matrix and RHS
% The product takes the form
% A[l,r x i,j] u[i,j] = F[l,r]
% l,r,i,j are the indices of the associated node in the element ei
AA   = zeros(G.NN);
FF   = zeros(G.NN,1);
idxb = zeros(G.NN,1);
for ei = 1 : length(EE)
  
  F  = zeros(G.N+1,G.N+1);
  A1 = zeros(G.N+1,G.N+1,G.N+1,G.N+1);
  A2 = zeros(G.N+1,G.N+1,G.N+1,G.N+1);
  A3 = zeros(G.N+1,G.N+1,G.N+1,G.N+1);
  A  = zeros(G.N+1,G.N+1,G.N+1,G.N+1);
  
  %#ok<*ALIGN>
  for l = 1:G.N+1, for r = 1:G.N+1     % RHS index
      for i = 1:G.N+1, for j = 1:G.N+1   % LHS index
          for m = 1:G.N+1, for n = 1:G.N+1 % GLL integration node index
              gp1 = [Ppkx2D(i,j,n,m) Ppky2D(i,j,n,m)]*EE(ei).Jinv(:,:,n,m);
              gp2 = [Ppkx2D(l,r,n,m) Ppky2D(l,r,n,m)]*EE(ei).Jinv(:,:,n,m);
              % Warning: A1 computation depends on mesh definition, this should be checked
              %         A1(l,r,i,j) = A1(l,r,i,j) + mu(EE(ei).XEK(l,r),EE(ei).YEK(l,r)).*EE(ei).AK(n,m).*(Ppkx2D(i,j,n,m).*Ppkx2D(l,r,n,m) + Ppky2D(i,j,n,m).*Ppky2D(l,r,n,m))/EE(ei).J(n,m);
              %         A1(l,r,i,j) = A1(l,r,i,j) + mu(EE(ei).XEK(l,r),EE(ei).YEK(l,r)).*EE(ei).AK(n,m).*(Px(i,j,n,m).*Px(l,r,n,m) + Py(i,j,n,m).*Py(l,r,n,m))/EE(ei).J(n,m);
              %         A1(l,r,i,j) = A1(l,r,i,j) + EE(ei).AK(n,m)*mu(EE(ei).XEK(n,m),EE(ei).YEK(n,m))*(gp1*gp2')/EE(ei).J(n,m);
              
              A1(l,r,i,j) = A1(l,r,i,j) + EE(ei).AK(n,m)*mu(EE(ei).XEK(n,m),EE(ei).YEK(n,m))*(gp1*gp2');
              % A2(l,r,i,j) = A2(l,r,i,j) + EE(ei).AK(n,m)*sigma(xk(m),xk(n)) * kronecker(i,l) * kronecker(j,r);
            end, end
          A2(l,r,l,r) = A2(l,r,l,r) + EE(ei).AK(l,r)*sigma(EE(ei).XEK(l,r),EE(ei).YEK(l,r)); % equivalent to commented line above
        end, end
      F(l,r) = EE(ei).AK(l,r)*f(EE(ei).XEK(l,r),EE(ei).YEK(l,r));
    end, end
  
  
  % BC  % CHECK WEIGHTS!
  for i = 1:G.N+1, for j = 1:G.N+1
      for m = 1:G.N+1
        %       %%% Top
        %       A3(G.N+1,m,i,j) = A3(G.N+1,m,i,j) + mu(EE(ei).XEK(G.N+1,m),EE(ei).YEK(G.N+1,m)).*EE(ei).dkt(m).*Ppky2D(i,j,G.N+1,m);
        %       %%% Bottom
        %       A3(1,m,i,j) = A3(1,m,i,j) + mu(EE(ei).XEK(1,m),EE(ei).YEK(1,m)).*EE(ei).dkb(m).*Ppky2D(i,j,1,m);
        %       %%% Right
        %       A3(m,G.N+1,i,j) = A3(m,G.N+1,i,j) + mu(EE(ei).XEK(m,G.N+1),EE(ei).YEK(m,G.N+1)).*EE(ei).dkr(m).*Ppkx2D(i,j,m,G.N+1);
        %       %%% Left
        %       A3(m,1,i,j) = A3(m,1,i,j) + mu(EE(ei).XEK(m,1),EE(ei).YEK(m,1)).*EE(ei).dkl(m).*Ppkx2D(i,j,m,1);

        
                %%% Top
                gpt = [Ppkx2D(i,j,G.N+1,m) Ppky2D(i,j,G.N+1,m)]*EE(ei).Jinv(:,:,G.N+1,m); % Jinv takes into account the fact that the actual derivative depends on the independent variable transformation
                A3(G.N+1,m,i,j) = A3(G.N+1,m,i,j) + EE(ei).dkb(m)*mu(EE(ei).XEK(G.N+1,m),EE(ei).YEK(G.N+1,m))*gpt(2); % The weights dk take into account the domain measurement modification
        
                %%% Bottom
                gpb = [Ppkx2D(i,j,1,m) Ppky2D(i,j,1,m)]*EE(ei).Jinv(:,:,1,m);
                A3(1,m,i,j) = A3(1,m,i,j) - EE(ei).dkt(m)*mu(EE(ei).XEK(1,m),EE(ei).YEK(1,m))*gpb(2);
        
                %%% Right
                gpr = [Ppkx2D(i,j,m,G.N+1) Ppky2D(i,j,m,G.N+1)]*EE(ei).Jinv(:,:,m,G.N+1);
                A3(m,G.N+1,i,j) = A3(m,G.N+1,i,j) + EE(ei).dkl(m)*mu(EE(ei).XEK(m,G.N+1),EE(ei).YEK(m,G.N+1))*gpr(1);
        
                %%% Left
                gpl = [Ppkx2D(i,j,m,1) Ppky2D(i,j,m,1)]*EE(ei).Jinv(:,:,m,1);
                A3(m,1,i,j) = A3(m,1,i,j) - EE(ei).dkr(m)*mu(EE(ei).XEK(m,1),EE(ei).YEK(m,1))*gpl(1); 
      end, end
  end
  
  % Sum contributes
  A = A1+A2-A3*0;
  
  
  % Recast 4-dimensional A into 2-dimensional AA and 2-dimensional F into
  % 1-dimensional FF
  for l = 1:G.N+1, for r = 1:G.N+1
      FF(EE(ei).gi(l,r)) = FF(EE(ei).gi(l,r)) + F(l,r);
      for i = 1:G.N+1, for j = 1:G.N+1
          AA(EE(ei).gi(l,r),EE(ei).gi(i,j)) = AA(EE(ei).gi(l,r),EE(ei).gi(i,j)) + A(l,r,i,j);
        end, end
      % Indexes of the boundary nodes
      idxb(EE(ei).gi(l,r),1) = EE(ei).bi(l,r);
      xn(EE(ei).gi(l,r),1)   = EE(ei).XEK(l,r);
      yn(EE(ei).gi(l,r),1)   = EE(ei).YEK(l,r);
    end, end
  
  
  
  % =======================================================================
  % Test different problem assembly
  for l = 1:G.N+1, for r = 1:G.N+1
      % linearize coordinates and nodes indexes (local, global, mortar flag)
      EE(ei).GI((l-1)*(G.N+1)+r,1) = EE(ei).gi(l,r);
      EE(ei).LI((l-1)*(G.N+1)+r,1) = EE(ei).li(l,r);
      EE(ei).MI((l-1)*(G.N+1)+r,1) = EE(ei).mi(l,r);
      EE(ei).BI((l-1)*(G.N+1)+r,1) = EE(ei).bi(l,r);
      EE(ei).XX((l-1)*(G.N+1)+r,1) = EE(ei).XEK(l,r);
      EE(ei).YY((l-1)*(G.N+1)+r,1) = EE(ei).YEK(l,r);
      % F and A
      EE(ei).FF((l-1)*(G.N+1)+r,1) = F(l,r);
      for i = 1:G.N+1, for j = 1:G.N+1
          EE(ei).AA((l-1)*(G.N+1)+r,(i-1)*(G.N+1)+j) = A(l,r,i,j);
        end, end
    end, end
  % Unify elements
  FFF(EE(ei).GI) = FFF(EE(ei).GI) + EE(ei).FF;
  AAA(EE(ei).GI,EE(ei).GI) = AAA(EE(ei).GI,EE(ei).GI) + EE(ei).AA;
  % =======================================================================
  
end



%% Assign boundary condition and solve system
% ub = g(xb,yb);
idxb = logical(idxb);
ub   = g(xn(idxb),yn(idxb));

gin = setdiff((1:G.NN),gb);
FF  = FF - AA(:,gb)*ub;
AA  = AA(gin,gin);
FF  = FF(gin);

Un      = zeros(G.NN,1);
Un(gin) = AA\FF;
Un(gb)  = ub;












%% Check against finite-differences Poisson solver
addpath poisson_test
figure
hold on
fh = testFdPoisson(G.xmin,G.xmax,f,g);
rmpath poisson_test
hold on
figure(fh)
plot3(xn,yn,Un,'.r','linewidth',2,'markersize',20)
plot3(xb,yb,ub,'ob','linewidth',2)
