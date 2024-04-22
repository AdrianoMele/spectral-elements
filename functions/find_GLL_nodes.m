function xk = find_GLL_nodes(N,x,Lp)
% xk = find_GLL_nodes(N) finds the N+1 GLL nodes of order N.
% xk = find_GLL_nodes(N,x,Lp) does the same, but does not compute x and
% L'(x) (where L(x) is the Legendre's polynomial of order N).

if nargin < 2, x = -1:1e-5:1; end
if nargin < 3, Lp = legendre_poly_derivative(N,x); end

% GLL quadrature nodes
xk = find_zeros(x,Lp); 
xk = (xk-flip(xk))/2; % symmetrize
xk = round(xk,10);
xk = [x(1) unique(xk) x(end)]; 