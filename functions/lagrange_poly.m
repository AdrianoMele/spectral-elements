function [L,xk] = lagrange_poly(N,x,xk,verbose)
% L = lagrange_poly(N,x,xk) computes Lagrange polynomials of order N over x
% with nodes xk. If the nodes are not specified, they are assumed to be
% equally spaced.

if nargin <= 2 || isempty(xk), xk = linspace(x(1),x(end),N+1); end
if nargin <= 3, verbose = 0; end

if (N ~= length(xk)-1), error('The number of nodes must be equal to N+1'); end

L = ones(N+1,length(x));
for j = 1 : N+1
  idx = setdiff(1:N+1,j);
  for m = idx
    L(j,:) = (x-xk(m))./(xk(j)-xk(m)) .* L(j,:);
  end
end

if verbose
  figure
  plot(x,L)
  hold on
  plot(xk,0*xk,'or'), grid
end

end

