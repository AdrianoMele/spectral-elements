%% Test gaussian integration
clearvars
close all
clc

%% Choose function for test
% f = @cos;
% f = @sin;
f = @(x) sin(x)+2*cos(x)+exp(x);
% f = @(x) 3*x.^4+x.^2+3*x;  % Polynomial integration should be exact for p < 2N-1

%% Initialize parameters
dx = 1e-3;
x  = -1:dx:1;
N  = 8; % Order of the integration

L = legendre_poly(N,x);
xk = find_zeros(x,L);
xk = unique(round(xk,4));
plot(x,L,xk,xk*0,'or')
pause

[L,~] = lagrange_poly(length(xk)-1,x,xk);
plot(x,L,xk,xk*0,'or');
pause

ak = sum(L,2)*dx;

a  = -4; b = 2;               % choose a different interval
tk = ((b-a)*xk + (a+b))/2;    % new nodes
dk = (b-a)*ak/2;              % new coefficients
x  = ((b-a)*x  + (a+b))/2;    % rescale x coordinate 
dx = x(2)-x(1);               % adapt x-step

%% Compute integral and show comparison
I = sum(f(x))*dx; % numerical
IG = f(tk)*dk;    % gaussian

% Plot interpolant
f_int = sum(L.*f(tk)');
plot(x,f(x),x,f_int,'--');
hold on
plot(tk,0*tk,'or')
title(['I = ' num2str(I) ' - IG = ' num2str(IG)])
grid

