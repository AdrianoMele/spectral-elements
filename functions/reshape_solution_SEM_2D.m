function [x,y,u] = reshape_solution_SEM_2D(xn,yn,un,G)
% [x,y,u] = reshape_solution_SEM_2D(xn,yn,un,G)

[xxn,ii] = sort(xn);
yyn      = yn(ii);
uun      = un(ii);

u = reshape(uun,G.MR*(G.N)+1,G.MC*(G.N)+1);
x = reshape(xxn,G.MR*(G.N)+1,G.MC*(G.N)+1);
y = reshape(yyn,G.MR*(G.N)+1,G.MC*(G.N)+1);
