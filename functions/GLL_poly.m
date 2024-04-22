function Psi = GLL_poly(N,x,xk,Lk)
% Psi = characteristic_GLL_poly(N,x,[xk],[Lk])
% Generate characteristic polynomials for GLL nodes xk evaluated at the
% points in x. Lk is the Legendre polynomial of order N evaluated at the
% nodes xk (optional).

if nargin < 3, xk = find_GLL_nodes(N); end
if nargin < 4, Lk = legendre_poly(N,xk); end
if min(x)<-1||max(x)>1, warning('The x vector should be in [-1,1]'); end

if iscolumn(x), x = x'; end

% Define new x to avoid numerical problems in the nodes
x1 = min([x xk])-1e-2:1e-3:max([x xk])+1e-2;
for i = 1:N+1
  x1 = x1(abs(x1-xk(i))>5e-2);
end
Lp = legendre_poly_derivative(N,x1);

% Compute polynomials
Psi = zeros(N+1,length(x1));
for j = 1 : N+1
  Psi(j,:) = -1/(N*(N+1)) * ( ((1-x1.^2).*Lp) ./ ((x1-xk(j)).*Lk(j)) );
end

% Manually add node values
x1 = [x1, xk];
Psi = [Psi, eye(N+1)];

% Sort everything
[x1,idx] = sort(x1);
Psi = Psi(:,idx);

% Go back to original x
Psi = interp1(x1,Psi',x,'spline')';












%%% Treat singular points
% % Psi(1,1)     = 1;
% % Psi(end,end) = 1;
% for j = 2:N
%   idx = find(x>=xk(j),1);
%   xtemp = [x(1:idx-nint)        xk(j)  x(idx+nint:end)];
%   ptemp = [Psi(j,1:idx-nint)    1      Psi(j,idx+nint:end)]; % Psi_i(xk(i)) = 1 by definition
%   Psi(j,:) = interp1(xtemp,ptemp,x,'spline');
% end