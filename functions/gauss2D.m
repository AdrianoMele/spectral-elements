function [ z ] = gauss2D(x,y,sigma,bound,varargin)
% Gaussian in 2D

if nargin > 4
  center = varargin{1};
else
  center = [0 0];
end

if length(x)~=length(y)
  error('x and y must be the same length.');
end

for i = 1:length(x)
  if(x(i)^2+y(i)^2) > bound^2
    z(i,1) = 0;
  else
    z(i,1) = 1/(2*pi*sigma) * exp(-((x(i)-center(1))^2+(y(i)-center(2))^2)/(2*sigma^2));
  end
end


end

