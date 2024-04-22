function z = find_zeros(x,f,verbose)
% z = find_zeros(x,f) finds the zeros of f(x) by linear interpolation of
% the samples.

if nargin < 3, verbose = 0; end
if iscolumn(x), x = x'; end
if length(x) ~= size(f,2), error('x and f must have the same number of samples'); end

z = [];
for j = 1 : size(f,1)
  y = f(j,:);
  i = find(diff(sign(y)));
  z = [z    x(i) - y(i).*(x(i+1)-x(i))./(y(i+1)-y(i))];
end

if verbose
  figure
  plot(x,f,z,0*z,'or')
  grid on
end
  