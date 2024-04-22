close all
addpath functions

figure, hold on, grid on

xr = [-1 +1 +1 -1];
yr = [-1 -1 +1 +1];
plot(xr([1:end 1]),yr([1:end 1]),'--sk','linewidth',2)


xe = [5 5   2   2.5];
ye = [1 0   1   2];
plot(xe([1:end 1]),ye([1:end 1]),'-sk','linewidth',2)

xk = find_GLL_nodes(3);
[XK,YK] = meshgrid(xk,xk);
[XEK,YEK,J,M] = affine_transform(xe,ye,XK,YK);

plot(XEK,YEK,'^r')
plot(XK,YK,'vr');
axis equal

rmpath functions