addpath functions
N = 5;
xk = find_GLL_nodes(N);
[XK,YK] = meshgrid(xk,xk);

idx = [];
xxk = zeros((N+1)^2,1); yyk = xxk;
for i = 1:N+1, for j = 1:N+1
%     idx = [idx; [(i-1)*(N+1)+j i j] ];
    xxk((i-1)*(N+1)+j) = XK(i,j);
    yyk((i-1)*(N+1)+j) = YK(i,j);
end, end

idxl = ((1:N+1)-1)*(N+1)+1;
idxr = ((1:N+1)-1)*(N+1)+(N+1);
idxb = (1-1)*(N+1) + (1:N+1);
idxt = (N+1-1)*(N+1) + (1:N+1);

idxbound = unique([idxt idxr idxb idxl]);
idxin    = setdiff(1:(N+1)^2,idxbound); 


close all
% plot(x(idxt),y(idxt),'ob');
xb = xxk(idxbound); yb = yyk(idxbound);
% plot(xb,yb,'ob');
ub = g(xb,yb);



rmpath functions