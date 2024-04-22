%% Test gaussian integration
clearvars
close all
clc

%% Choose function for test
% f = @cos;
% f = @sin;
f = @(x) sin(x)+2*cos(x)+exp(x);
% f = @(x) x.^2+3*x;  % Polynomial integration should be exact for p < 2N-1

%% Initialize parameters
dx = 1e-2;
a = -1; b = 1;
x  = a:dx:b;
N  = 4; % Order of the integration

xk = linspace(-1,1,N+1);
xk = ((b-a)*xk + (a+b))/2;

[L,~] = lagrange_poly(length(xk)-1,x,xk);
ak = sum(L,2)*dx;
% dk = (b-a)*ak/2; % automatically verified

%% Compute integral and show comparison
I = sum(f(x))*dx; % numerical
IG = f(xk)*ak;    % gaussian

% Plot interpolant
f_int = sum(L.*f(xk)');
plot(x,f(x),x,f_int,'--');
hold on
plot(xk,0*xk,'or')
title(['I = ' num2str(I) ' - IG = ' num2str(IG)])
grid

