function fh = testFdPoisson(a,b,f,g)

% f = inline('-5*pi^2*sin(pi*x).*cos(2*pi*y)');
% g = inline('sin(pi*x).*cos(2*pi*y)');

% f = inline('-(x.^2+y.^2)');
% g = inline('0*x.*y');

ff = @(x,y) -f(x,y);
[u,x,y] = fd2poisson(ff,g,a,b,39);
% u = reshape(u,length(y),length(x));
h = x(1,2) - x(1,1);

%% Plot solution
% fh = figure; 
set(gcf,'DefaultAxesFontSize',8,'PaperPosition', [0 0 3.5 3.5]);
mesh(x,y,u), colormap([0 0 0]), xlabel('x'), ylabel('y'), zlabel('u(x,y)');
title(strcat('Numerical Solution to Poisson Equation, h=',num2str(h)));
fh = gcf;




return
%% Plot error
figure, set(gcf,'DefaultAxesFontSize',8,'PaperPosition', [0 0 3.5 3.5]),
mesh(x,y,u-g(x,y)), colormap([0 0 0]), xlabel('x'), ylabel('y'), zlabel('Error'),
title(strcat('Error, h=',num2str(h)));