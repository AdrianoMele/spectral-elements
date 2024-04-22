function [AA,Un] = SEM_2D_operator(G,EE,g,mu,sigma,Ppkx2D,Ppky2D,xb,yb,gb,U0)
% [AA,[Un]] = SEM_2D_operator(G,EE,g,mu,sigma,xb,yb,gb,[Ppkx2D],[Ppky2D],[U0],[wb])
% 
% Computes the SEM_2D operator for a generic linear differential operator 
% of order 2 defined as
%   Lu = -div(mu(x,y)*grad(u(x,y))) + sigma(x,y)*u
% 
% u >= u0 > 0, s >= 0 are needed for the bilinear form to be coercive
% (See sample problem in Quarteroni et al. p.245)
% 
% The operator is returned as a matrix AA, associated to the nodes
% specified in EE.
% 
% G is the geometry definition and EE the mesh struct 
% (see also initmeshSEM)
% 
% g, f, mu and sigma are function handles (defined above)
% 
% Ppkx2D,Ppky2D are the derivatives of the GLL polynomials (optional, see
% GLL_poly_2D_gradient)
% 
% xb,yb,gb are the coordinates and the global indexes of the boundary nodes
% 
% U0 is the vector to which the operator must be applied (optional). If it
% is specified, the operator is applied to it and the value at the nodes 
% specified in EE is returned in Un (optional output argument).

wb = waitbar(0,'Building discrete operator...');

if exist('Ppkx2D','var')==0
  [xk,~,Lk,~,~,~] = init_ref_element(G.N);
  [Ppkx2D,Ppky2D] = GLL_poly_2D_gradient(G.N,xk,xk,xk,Lk);
end

% Stiffness matrix and RHS
% The product takes the form A[l,r x i,j] u[i,j] = F[l,r]
% l,r and i,j are the indices of the associated node in the element #ei
AA   = zeros(G.NN);
for ei = 1 : length(EE)
  
  A1 = zeros(G.N+1,G.N+1,G.N+1,G.N+1);
  A2 = zeros(G.N+1,G.N+1,G.N+1,G.N+1);
  A3 = zeros(G.N+1,G.N+1,G.N+1,G.N+1);
  
  %#ok<*ALIGN>
  for l = 1:G.N+1, for r = 1:G.N+1     % RHS index
    for i = 1:G.N+1, for j = 1:G.N+1   % LHS index
      for m = 1:G.N+1, for n = 1:G.N+1 % GLL integration node index
        gp1 = [Ppkx2D(i,j,n,m) Ppky2D(i,j,n,m)]*EE(ei).Jinv(:,:,n,m);
        gp2 = [Ppkx2D(l,r,n,m) Ppky2D(l,r,n,m)]*EE(ei).Jinv(:,:,n,m);
        % Warning: A1 computation depends on mesh definition, this should be checked          
        A1(l,r,i,j) = A1(l,r,i,j) + EE(ei).AK(n,m)*mu(EE(ei).XEK(n,m),EE(ei).YEK(n,m))*(gp1*gp2');
%         A2(l,r,i,j) = A2(l,r,i,j) + EE(ei).AK(n,m)*sigma(xk(m),xk(n)) * kronecker(i,l) * kronecker(j,r); % equivalent formulation for A2 (see below)
      end, end
      A2(l,r,l,r) = A2(l,r,l,r) + EE(ei).AK(l,r)*sigma(EE(ei).XEK(l,r),EE(ei).YEK(l,r)); % equivalent to commented line above
    end, end
  end, end
 

%   Treat boundary:  
%   This part is needed ONLY at the domain boundary. Since we are using
%   Dirichlet BC, it can be safely ignored. In the inner points, the 
%   condition on one element is compensated by an opposite condition on 
%   the other. In the future, this part could be useful to treat the 
%   general BC case.
% 
%   for i = 1:G.N+1, for j = 1:G.N+1
%     for m = 1:G.N+1
%       %%% Top
%       gpt = [Ppkx2D(i,j,G.N+1,m) Ppky2D(i,j,G.N+1,m)]*EE(ei).Jinv(:,:,G.N+1,m); % Jinv takes into account the fact that the actual derivative depends on the independent variable transformation
%       A3(G.N+1,m,i,j) = A3(G.N+1,m,i,j) + EE(ei).dkb(m)*mu(EE(ei).XEK(G.N+1,m),EE(ei).YEK(G.N+1,m))*gpt(2); % The weights dk take into account the domain measurement modification       
%       %%% Bottom
%       gpb = [Ppkx2D(i,j,1,m) Ppky2D(i,j,1,m)]*EE(ei).Jinv(:,:,1,m);
%       A3(1,m,i,j) = A3(1,m,i,j) - EE(ei).dkt(m)*mu(EE(ei).XEK(1,m),EE(ei).YEK(1,m))*gpb(2);      
%       %%% Right
%       gpr = [Ppkx2D(i,j,m,G.N+1) Ppky2D(i,j,m,G.N+1)]*EE(ei).Jinv(:,:,m,G.N+1);
%       A3(m,G.N+1,i,j) = A3(m,G.N+1,i,j) + EE(ei).dkl(m)*mu(EE(ei).XEK(m,G.N+1),EE(ei).YEK(m,G.N+1))*gpr(1);        
%       %%% Left
%       gpl = [Ppkx2D(i,j,m,1) Ppky2D(i,j,m,1)]*EE(ei).Jinv(:,:,m,1);
%       A3(m,1,i,j) = A3(m,1,i,j) - EE(ei).dkr(m)*mu(EE(ei).XEK(m,1),EE(ei).YEK(m,1))*gpl(1); 
%     end, end
%   end
 
  % Recast 4D A into 2D AA
  A = A1+A2-A3;
  for l = 1:G.N+1, for r = 1:G.N+1
    for i = 1:G.N+1, for j = 1:G.N+1
      AA(EE(ei).gi(l,r),EE(ei).gi(i,j)) = AA(EE(ei).gi(l,r),EE(ei).gi(i,j)) + A(l,r,i,j);
    end, end
  end, end  

  waitbar(ei/numel(EE),wb)
end

close(wb)
%% Apply operator to assigned U if required
if nargout>1 && exist('U0','var')~=0
  gin = setdiff((1:G.NN),gb);
  Un(gin) = AA(gin,gin)*U0(gin);
  Un(gb)  = g(xb,yb);
end