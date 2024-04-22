function [xk,ak,Lk,XK,YK,AK] = init_ref_element(N)
% [xk,ak,XK,YK,LK,AK] = init_ref_element(N)
% Initializes reference element for GLL integration.
% The 1D reference element is [-1,1]. xk are the GLL nodes, LK are the 
% Legendre Polynomials evaluated at the nodes and ak the numerical
% integration coefficients.
% The 2D reference element is [-1,1]x[-1,1]. XK and YK are the coordinates 
% of the GLL nodes, AK are the GLL quadrature integration coefficients. 

% GLL quadrature nodes
xk = find_GLL_nodes(N);
[XK,YK] = meshgrid(xk,xk);

% Integration coefficients
Lk = legendre_poly(N,xk);
ak = 2/(N*(N+1)) ./ (Lk'.^2);   % Analytical formula
AK = ak*ak';