function c = binomial(n,k)
% Returns binomial coefficient n over k.

% if ~isinteger(n) || ~isinteger(k), error('Input values must be integers'); end
if k > n, error('k should be less or equal to n'); end

c = factorial(n)/(factorial(k)*factorial(n-k));