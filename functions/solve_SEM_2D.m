function Un = solve_SEM_2D(G,EE,g,f,mu,sigma,xb,yb,gb,AA)
% Un = solve_SEM_2D(G,EE,g,f,mu,sigma,xb,yb,gb,[AA],[EE])
% 
% Assembles and solves the 2D SEM system for the problem
%           Lu(x,y) = f(x,y) in [a,b]x[c,d]
% with BC   u = g(x,y) on the boundary
% 
% L is a generic linear differential operator of order 2 defined as
%   Lu = -div(mu(x,y)*grad(u(x,y))) + sigma(x,y)*u
% 
% u >= u0 > 0, s >= 0 are needed for the bilinear form to be coercive
% (See sample problem in Quarteroni et al. p.245)
% 
% G is the geometry definition and EE the mesh struct 
% (see also initmeshSEM)
% 
% g, f, mu and sigma are function handles (defined above)
% 
% xb,yb,gb are the coordinates and the global indexes of the boundary nodes
% 
% The optional parameter AA can be use to recover a matrix 
% operator built in a previous run.
% 
% Un are the values of the solution at the nodes specified in EE.

if ~exist('AA','var'), AA = SEM_2D_operator(G,EE,g,mu,sigma); end
FF = SEM_2D_RHS(G,EE,f);

Un     = zeros(G.NN,1);
Un(gb) = g(xb,yb);

gin = setdiff((1:G.NN),gb);
FF  = FF - AA(:,gb)*Un(gb);
Un(gin) = AA(gin,gin)\FF(gin);
