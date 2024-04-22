function Psi_p = GLL_poly_derivatives(N,x,xk,Lk,eps,verbose)
% Psi_p = GLL_poly_derivatives(N,x,[xk],[Lk],[eps=1e-4],[verbose=0])
% Generate derivatives of the characteristic polynomials of order N for GLL
% nodes. The derivatives are computed numerically.

if nargin < 3, xk = find_GLL_nodes(N); end
if nargin < 4 || isempty(Lk),  Lk = legendre_poly(N,xk); end
if nargin < 5 || isempty(eps), eps = 1e-4; end
if nargin < 6, verbose = 0; end

warning off
Psi1 = GLL_poly(N,x-eps,xk,Lk);
Psi2 = GLL_poly(N,x+eps,xk,Lk);
warning on

Psi_p = (Psi2-Psi1)/(2*eps);

if verbose
  figure
  for i = 1 : N+1
    plot(x,Psi1(i,:),'k', x(1:end-1),diff(Psi1(i,:))./diff(x),'b', x,Psi_p(i,:),'.r')
    grid on
    pause
  end
end