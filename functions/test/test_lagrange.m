function test_lagrange(N,x)
% Test lagrange polynomials derivatives computation

close all
clc

if nargin<2
  x = -1:1e-4:1;
end

P  = lagrange_poly(N,x);                      % ok
P2 = lagrange_poly_derivatives(N,x);          % ok
P1 = cumsum(P2*1e-4,2) + repmat(P(:,1),1,[]); % la primitiva è definita a meno di una costante

figure
for i = 1 : N+1
  subplot(121) % integrale numerico vs polinomio
  plot(x,P(i,:),x,P1(i,:),'--')
  grid on
  
  subplot(122) % derivata numerica vs analitica
  plot(x,P2(i,:),x(1:end-1),diff(P(i,:))./diff(x), '--')
  grid on
  
  pause
end