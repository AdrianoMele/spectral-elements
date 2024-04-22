function [L,nk] = lagrange_poly_derivatives(N,x,xk,verbose)
% L = lagrange_poly(N,x,xk) computes the derivatives wrt x of the
% legendre polynomials of order N over x with nodes xk. If the nodes are
% not specified, they are assumed to be equally spaced.

if nargin <= 2 || isempty(xk), xk = linspace(x(1),x(end),N+1); end
if nargin <= 3, verbose = 0; end

L = zeros(N+1,length(x));
for j = 1 : N+1
  idx1 = setdiff(1:N+1,j);
  for n = idx1
    idx2 = setdiff(idx1,n);
    l = ones(1,length(x));
    for m = idx2
      l = (x-xk(m))./(xk(j)-xk(m)) .* l;
    end
    L(j,:) = L(j,:) + l./(xk(j)-xk(n));
  end
end

% Find nodes
nk = [x(1) find_zeros(x,L) x(end)];
nk = round(nk,5);
nk = unique(sort(nk));

if verbose
  figure
  plot(x,L)
  hold on
  plot(nk,0*nk,'or'), grid
end
end


