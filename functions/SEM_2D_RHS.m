function FF = SEM_2D_RHS(G,EE,f)
% FF = SEM_2D_RHS(G,EE,f)
%
% Computes RHS for problem in solve_SEM_2D

FF   = zeros(G.NN,1);
for ei = 1 : length(EE)
  F  = zeros(G.N+1,G.N+1);
  %#ok<*ALIGN>
  for l = 1:G.N+1, for r = 1:G.N+1     % RHS index
    F(l,r) = EE(ei).AK(l,r)*f(EE(ei).XEK(l,r),EE(ei).YEK(l,r));
  end, end
  for l = 1:G.N+1, for r = 1:G.N+1
    FF(EE(ei).gi(l,r)) = FF(EE(ei).gi(l,r)) + F(l,r);    
  end, end
  
end