%% Test gauss-lobatto-legendre integration in 2D
addpath ./functions
clearvars
close all
clc

verbosity = 1;


%% Choose function for test
% f = @cos;
f = @(x,y) cos(x).*cos(y);
% f = @(x,y) sin(x)+2*cos(y)+exp(-x+y.^2);
% f = @(x,y) x.^2-3*x+y.^4-y.^3+x.*y;  % Polynomial integration should be exact for p < 2N-1


%% Initialize parameters
dx = 5e-2;    dy = dx*2;
x  = -1:dx:1; y = -1:dy:1;
N  = 5; % Order of the integration

% GLL quadrature nodes
xk = find_GLL_nodes(N);

% Integration coefficients
Lk = legendre_poly(N,xk);
ak = 2/(N*(N+1)) ./ (Lk'.^2);   % Analytical formula
ajk = ak*ak';


%% Find interpolant
% Characteristic polynomials for GLL nodes
% Psi   = characteristic_GLL_poly(N,x,xk,Lk);
% Psi2D = zeros(N+1,N+1,length(y),length(x));
% for j = 1 : N+1
%   for k = 1 : N+1
%     Psi2D(j,k,:,:) = Psi(j,:)'*Psi(k,:);
%   end
% end

Psi2D = GLL_poly_2D(N,x,y,xk,Lk);

if verbosity
  [xx,yy] = meshgrid(x,y);
  figure(1)
  for j = 1:N+1
    for k = 1:N+1
      subplot(N+1,N+1,((j-1)*(N+1))+k)
      surf(xx,yy,squeeze(Psi2D(j,k,:,:)));
    end
  end
end

%% Change domain definition

% X-direction
a   = -1; b = 4;               % choose a different interval
txk = ((b-a)*xk + (a+b))/2;    % new nodes
dxk = (b-a)*ak/2;              % new coefficients
x   = ((b-a)*x  + (a+b))/2;    % rescale x coordinate
dx  = x(2)-x(1);               % adapt x-step
% Y-direction
c   = -1; d = 3;               % choose a different interval
tyk = ((d-c)*xk + (c+d))/2;    % new nodes
dyk = (d-c)*ak/2;              % new coefficients
y   = ((d-c)*y  + (c+d))/2;    % rescale y coordinate
dy  = y(2)-y(1);               % adapt y-step

djk = dxk*dyk';                % new weights


%% 2D function
% Function on grid points and numerical integral
[xx,yy] = meshgrid(x,y);
F = xx*0;
for c = 1 : length(x)
  for r = 1 : length(y)
    F(r,c) = f(xx(r,c),yy(r,c));
  end
end
I1 = integral2(f, x(1),x(end),y(1),y(end));

% Function on nodes and GLL integral
[xxk,yyk] = meshgrid(txk,tyk);
Fjk = zeros(N+1,N+1);
for c = 1 : length(txk)
  for r = 1 : length(tyk)
    Fjk(r,c) = f(xxk(r,c),yyk(r,c));
  end
end
I2 = sum(sum(djk.*Fjk));

if verbosity
  figure(2)
  plot3(xx,yy,F,'ob');
  hold on, grid on
  plot3(xxk,yyk,xxk*0,'+k','linewidth', 2);
  title(['I1 = ' num2str(I1) ' - I2 = ' num2str(I2)])
else
  disp(['I1 = ' num2str(I1) ' - I2 = ' num2str(I2)]);
end

%% Plot interpolant
if verbosity
  Fint = xx*0;
  for j = 1:N+1
    for k = 1:N+1
      Fint = Fint + squeeze(Fjk(j,k).*Psi2D(j,k,:,:));
    end
  end
  figure(2)
  plot3(xx,yy,Fint,'.r');
end





