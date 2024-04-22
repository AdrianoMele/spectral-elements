function test_GLL_poly_2D(N,x,y,xk,Lk)
% test_GLL_poly_2D(N,[x],[y],[xk],[Lk])
% Test the computation of the GLL characteristic polynomials of order N and 
% of their gradient over the grid defined by meshgrid(x,y).

if nargin < 2, x  = linspace(-1,1,100);  end
if nargin < 3, y  = linspace(-1,1,100);  end
if nargin < 4, xk = find_GLL_nodes(N);   end
if nargin < 5, Lk = legendre_poly(N,xk); end

Psi2D = GLL_poly_2D(N,x,y,xk,Lk);
[Px2D,Py2D,XG,YG,XK,YK] = GLL_poly_2D_gradient(N,x,y,xk,Lk);

% Check GLL_Poly
for i = 1 : length(xk)
  for j = 1 : length(xk)
    subplot(1,2,1)
    mesh(XG,YG,squeeze(Psi2D(i,j,:,:)))
    hold on
    plot3(XK(i,j),YK(i,j),1,'or')
    hold off
    xlabel('x'), ylabel('y'),zlabel('\Psi(x,y)')
%     subplot(1,3,2)
%     mesh(XG,YG,squeeze(Ppkx2D(i,j,:,:)))
%     xlabel('x'), ylabel('y'),zlabel('\partial\Psi / \partial x')
%     subplot(1,3,3)
%     mesh(XG,YG,squeeze(Ppky2D(i,j,:,:)))
%     xlabel('x'), ylabel('y'),zlabel('\partial\Psi / \partial y')
    
    subplot(1,2,2)
    contour(XG,YG,squeeze(Psi2D(i,j,:,:)),100);
    hold on, axis equal
    quiver(XG,YG,squeeze(Px2D(i,j,:,:)),squeeze(Py2D(i,j,:,:)));
    plot(XK,YK,'xr')
    hold off
    xlabel('x'), ylabel('y')
    title('\nabla\Psi(x,y) vector field and \Psi contour lines')
        
    shg
    pause
  end
end


dx = x(2)-x(1);
dy = y(2)-y(1);
for i = 1:N
  for j = 1:N
    [ppx,ppy]=gradient(squeeze(Psi2D(i,j,:,:)),dx,dy);
    
    subplot(121)
    plot3(XG,YG,squeeze(Px2D(i,j,:,:)),'ob')
    hold on
    plot3(XG,YG,ppx,'.r')
    hold off
    title('\Psi_x')
    
    
    subplot(122)
    plot3(XG,YG,squeeze(Py2D(i,j,:,:)),'ob')
    hold on
    plot3(XG,YG,ppy,'.r') 
    hold off
    title('\Psi_y')
    
    pause
  end
end