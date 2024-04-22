%% Test gauss-lobatto-legendre integration
clearvars
close all
clc


%% Choose function for test
% f = @cos;
% f = @sin;
f = @(x) sin(x)+2*cos(x)+exp(-x+2);
% f = @(x) x.^2+3*x;  % Polynomial integration should be exact for p < 2N-1


%% Initialize parameters
dx = 1e-3;
x  = -1:dx:1;
N  = 8; % Order of the integration

% Legendre polynomials
L  = legendre_poly(N,x);

% GLL quadrature nodes
xk = find_GLL_nodes(N); 

% Characteristic polynomials for GLL nodes
Lk = legendre_poly(N,xk);
Psi = GLL_poly(N,x,xk,Lk);
% plot(x,Psi,'-^',xk,xk*0,'or'), grid on, hold on % For N=4, check Quarteroni, p. 233


%% Compute integration coefficients
ak = sum(Psi,2)*dx;           % Numerical
% ak = 2/(N*(N+1)) ./ (Lk'.^2);   % Analytical


%% Change domain definition
a  = 1; b = 10;               % choose a different interval
tk = ((b-a)*xk + (a+b))/2;    % new nodes
dk = (b-a)*ak/2;              % new coefficients
x  = ((b-a)*x  + (a+b))/2;    % rescale x coordinate 
dx = x(2)-x(1);               % adapt x-step


%% Compute integral and show comparison
I = sum(f(x))*dx; % numerical
IG = f(tk)*dk;    % gaussian


%% Plot interpolant and results
f_int = sum(Psi.*f(tk)');
plot(x,f(x),x,f_int,'--');
hold on
plot(tk,0*tk,'or')
title(['I = ' num2str(I) ' - IG = ' num2str(IG)])
grid

