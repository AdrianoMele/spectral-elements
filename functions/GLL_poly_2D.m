function [Psi2D,X,Y,XK,YK] = GLL_poly_2D(N,x,y,xk,Lk)
% Psi = GLL_poly_2D(N,x,y,[xk],[Lk])
% [Psi,X,Y] = GLL_poly_2D(N,x,y,[xk],[Lk])
% [Psi,X,Y,XK,YK] = GLL_poly_2D(N,x,y,[xk],[Lk])
% Generate 2D characteristic polynomials for GLL nodes over the grid 
% [X,Y] = meshgrid(x,y) for the nodes [XK,YK] = meshgrid(xk). Lk is the 1D 
% Legendre polynomial of order N evaluated at the nodes xk (optional).
%
% Psi is a [N+1]x[N+1]x[Ny]x[Nx] matrix, X and Y are [Ny]x[Nx] matrices and
% XK and YK are [N+1]x[N+1] matrices. 
% The value of Psi(n,m,:,:) is 1 at the node (XK(n,m),YK(n,m)) and zero at 
% the others.
%
% See also GLL_poly

if nargin < 4, xk = find_GLL_nodes(N); end
if nargin < 5, Lk = legendre_poly(N,xk); end

if iscolumn(x), x = x'; end
if iscolumn(x), y = y'; end

Psix  = GLL_poly(N,x,xk,Lk);
Psiy  = GLL_poly(N,y,xk,Lk);

Psi2D = zeros(N+1,N+1,length(y),length(x));
for j = 1 : N+1
  for k = 1 : N+1
    Psi2D(j,k,:,:) = Psiy(j,:)'*Psix(k,:);
  end
end

if nargout > 1, [X,Y]   = meshgrid(x,y);   end
if nargout > 3, [XK,YK] = meshgrid(xk,xk); end

end

