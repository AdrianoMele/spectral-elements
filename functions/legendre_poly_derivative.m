function L = legendre_poly_derivative(N,x)
% L = legendre_poly(N,x) returns the derivative of the  Legendre polynomial 
% of order N over the vector x.

if iscolumn(x), x = x'; end


L = zeros(size(x));
% k = 0
L = L + binomial(N,0)^2* (N*(x-1).^(N-1));
% k = 1, ..., N-1
for k = 1 : N-1
  L = L + binomial(N,k)^2* ((N-k)*((x-1).^(N-k-1)).*((x+1).^k) + (k*(x-1).^(N-k)).*((x+1).^(k-1)));
end
% k = N
L = L + binomial(N,N)^2* (N.*((x+1).^(N-1)));

L = 1/(2^N) * L;



% switch N
%   case 0
%     L = ones(size(x));
%   case 1
%     L = x;
%   case 2
%     L = (3*x.^2-1)/2;
%   case 3
%     L = (5*x.^3 - 3*x)/2;
%   case 4
%     L = (35*x.^4 - 30*x.^2 + 3)/8;
%   case 5 
%     L = (63*x.^5 - 70*x.^3 + 15*x)/8;
%   case 6
%     L = (231*x.^6 - 315*x.^4 + 105*x.^2 - 5)/16;
%   case 7
%     L = (429*x.^7 - 693*x.^5 + 315*x.^3 - 35*x)/16;
%   case 8
%     L = (6435*x.^8 - 12012*x.^6 + 6930*x.^4 - 1260*x.^2 + 35)/128;
%   case 9
%     L = (12155*x.^9 - 25740*x.^7 + 18018*x.^5 - 4620*x.^3 + 315*x)/128;
%   case 10
%     L = (46189*x.^10 - 109395*x.^8 + 90090*x.^6 - 30030*x.^4 + 3465*x.^2 - 63)/256;
%   otherwise 
%     error('Only the first 10 Legendre polynomials are available');
% end


