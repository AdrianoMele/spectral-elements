function [Px,Py,X,Y,XK,YK] = GLL_poly_2D_gradient(N,x,y,xk,Lk)
% [Px,Py,X,Y,XK,YK] = GLL_poly_2D_gradient(N,x,y,[xk],[Lk])
% Computes the gradient of the 2D characteristic polynomials for GLL
% integration. The gradient is evaluated at the points specified by the x
% and y vectors. Lk is the Legendre polynomial of order N evaluated at the
% nodes xk (optional)
% 
% Px and Py are [N+1]x[N+1]x[Ny]x[Nx] matrices, X and Y are [Ny]x[Nx] 
% matrices and XK and YK are [N+1]x[N+1] matrices. 
% The value of Psi(n,m,:,:) is 1 at the node (XK(n,m),YK(n,m)) and zero at 
% the others (i.e. Psi(n,m,:,:) is the node function associated to 
% (XK(n,m),YK(n,m)).
%
% See also GLL_poly, GLL_poly_derivatives


if nargin < 4, xk = find_GLL_nodes(N);   end
if nargin < 5, Lk = legendre_poly(N,xk); end

if iscolumn(x), x = x'; end
if iscolumn(y), y = y'; end

nx = length(x);
ny = length(y);

Psix   = GLL_poly(N,x,xk,Lk);
Psiy   = GLL_poly(N,y,xk,Lk);
Psi_px = GLL_poly_derivatives(N,x,xk,Lk);
Psi_py = GLL_poly_derivatives(N,y,xk,Lk);

Px = zeros(N+1,N+1,ny,nx);
Py = zeros(N+1,N+1,ny,nx);
for j = 1:N+1
  for k = 1:N+1
    Px(j,k,:,:) = Psiy(j,:)'*Psi_px(k,:);
    Py(j,k,:,:) = Psi_py(j,:)'*Psix(k,:);
  end
end

if nargout>2, [X,Y]   = meshgrid(x,y); end
if nargout>4, [XK,YK] = meshgrid(xk,xk); end