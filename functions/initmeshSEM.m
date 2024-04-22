function [EE,xn,yn,gi,xb,yb,gb,xm,ym,gm] = initmeshSEM(G,verbose)
% [EE,xn,yn,gi,bi,mi,xb,yb,gb,xm,ym,gm] = initmeshSEM(G,[verbose=0])
%
% G is the mesh geometry definition, containing:
%     G.N     = GLL integration order
%     G.MR    = no. of mesh rows
%     G.MC    = no. of mesh columns
%     G.xmin  = grid minimum x
%     G.xmax  = grid maximum x
%     G.ymin  = grid minimum y
%     G.ymax  = grid maximum y
%
% EE is a structure containing for each element the following fields:
%     EE(ei).xve  = vertices x-coordinates;
%     EE(ei).yve  = vertices y-coordinates;
%     EE(ei).M    = element measure;
%     EE(ei).XEK  = GLL integration nodes, x coordinates;
%     EE(ei).YEK  = GLL integration nodes, y coordinates;
%     EE(ei).J    = GLL integration nodes, Jacobian;
%     EE(ei).AK   = GLL integration weights;
%     EE(ei).li   = nodes local indexes;
%     EE(ei).gi   = nodes global indexes;
%     EE(ei).bi   = nodes boundary flags;
%     EE(ei).Jinv = GLL integration nodes, full inverse jacobian matrix;
%     EE(ei).dkt  = GLL line integration weights, top edge;
%     EE(ei).dkb  = GLL line integration weights, bottom edge;
%     EE(ei).dkr  = GLL line integration weights, left edge;
%     EE(ei).dkl  = GLL line integration weights, right edge;
%     EE(ei).mi   = 1 if the node is a mortar node (master), -1 if slave
% 
% xn,yn are the nodes coordinates
% gi are the global indices associated to the nodes
% xb,yb,gb are the coordinates and indexes of the boundary nodes
% xm,ym,gm are the coordinates and indexes of the mortar nodes
%  
% See also affine_transform.m

%%% Define a structure to order the elements. Informations needed:
% Index of the element
% Vertexes coordinates
% Area of the element
% Local indexes of the nodes
% Global indexes of the nodes
% Nodes coordinates
% Integration weights
% Jacobian of the affine transformation

% The global numbering of the elements follows this convention:
%
%     NELX  2*NELY   ... NELX*NELY
%      ...   ...     ...   ... 
% ^     2    ...     ...   ... 
% |     1   NELY+1   ...   ... 
% --->

% The local numbering of GLL nodes follows this convention:
% 
%  	(1,NGLL)...	(NGLL,NGLL)
% ^   	...	...	... 
% |	(1,1)	...	(NGLL,1)
% --->


if nargin==1, verbose = 0; end
if verbose, fmesh = figure; end

% Initialize reference element
[~,ak,~,XK,YK,AK] = init_ref_element(G.N);

% Define coordinates of the grid nodes.
% Nodes are ordered from bottom-left to top-right.
xm = linspace(G.xmin,G.xmax,G.MC+1);
ym = linspace(G.ymin,G.ymax,G.MR+1);
[XM,YM] = meshgrid(xm,ym);

% Initialize global element index temporary matrix
gi =  zeros(G.N+1);
for i = 1 : G.N+1
  for j = 1 : G.N+1
    gi(i,j) = (i-1)*(G.N+1)+j;
  end
end
gi_curr = 1;

% Build elements information structure
EE = [];
for mc = 1:G.MC
  for mr = 1:G.MR
    
    % Current element index
    ei = (mc-1)*G.MR+mr;
    
    % Vertexes of the element
    a = XM(mr,mc);
    b = XM(mr,mc+1);
    c = YM(mr,mc);
    d = YM(mr+1,mc);
    
    % Nodes
    xve = [a b b a]; yve = [c c d d];
    [XEK,YEK,J,M,JJ] = affine_transform(xve,yve,XK,YK);
    
    % Local index of the nodes
    li = zeros(G.N+1,G.N+1);
    for i = 1 : G.N+1
      for j = 1 : G.N+1
        li(i,j) = (i-1)*(G.N+1)+j;
      end
    end
    
    % Global index of the nodes
    gi = zeros(G.N+1);
    if mc == 1
      idxc = 1:G.N+1;
    else
      idxc = 2:G.N+1;
      gi(1:G.N+1,1) = EE(ei-G.MR).gi(1:G.N+1,G.N+1);
    end
    if mr == 1
      idxr = 1:G.N+1;
    else
      idxr = 2:G.N+1;
      gi(1,1:G.N+1) = EE(ei-1).gi(G.N+1,1:G.N+1);
    end
    
    for i = idxr
      for j = idxc
        gi(i,j) = gi_curr; gi_curr = gi_curr+1;
      end
    end

    % Identify mortar nodes
    mi = false(G.N+1);
    mi(:,1)   = true;
    mi(1,:)   = true;
    mi(:,end) = true;
    mi(end,:) = true;
        
    % Identify boundary nodes
    bi = false(G.N+1);
    if mc == 1
      bi(1:G.N+1,1) = true;
      mi(1:G.N+1,1) = false;
    end
    if mc == G.MC
      bi(1:G.N+1,G.N+1) = true;
      mi(1:G.N+1,G.N+1) = false;
    end
    
    if mr == 1
      bi(1,1:G.N+1) = true;
      mi(1,1:G.N+1) = false;
    end
    if mr == G.MR
      bi(G.N+1,1:G.N+1) = true;
      mi(G.N+1,1:G.N+1) = false;
    end
        
    % Store informations for each element in a struct
    EE(ei).xve  = xve;
    EE(ei).yve  = yve;
    EE(ei).M    = M;
    EE(ei).XEK  = XEK;
    EE(ei).YEK  = YEK;
    EE(ei).J    = J;
%     EE(ei).AK   = AK*M/4;
    EE(ei).AK   = AK.*EE(ei).J; % equivalent to the line above
    EE(ei).li   = li;
    EE(ei).gi   = gi;
    EE(ei).bi   = bi;
    EE(ei).Jinv = JJ;
    EE(ei).dkt  = (XEK(1,end)-XEK(1,1))/2*ak;
    EE(ei).dkb  = (XEK(end,end)-XEK(end,1))/2*ak;
    EE(ei).dkl  = (YEK(end,1)-YEK(1,1))/2*ak;
    EE(ei).dkr  = (YEK(end,end)-YEK(1,end))/2*ak;
    EE(ei).mi   = mi;
    
    % Plot current element
    if verbose
      figure(fmesh);
      hold on, grid on
      plot(xve([1:4 1]), yve([1:4 1]),'-ob','linewidth',1)
      axis equal
      for i = 1:G.N+1
        for j = 1:G.N+1
          text(XEK(i,j),YEK(i,j),int2str(gi(i,j)),'fontsize',20)
        end
      end
      plot(XEK,YEK,'xk')
%       title(['el ' int2str(ei) ' - row ' int2str(mr) ' - col ' int2str(mc)])
%       pause
    end
    
  end
end

gi = [];
bi = [];
mi = [];
xn = [];
yn = [];

for ei = 1 : length(EE)
  gi = [gi; EE(ei).gi(:)];
  bi = [bi; logical(EE(ei).bi(:))];
  mi = [mi; logical(EE(ei).mi(:))];
  xn = [xn; EE(ei).XEK(:)];
  yn  = [yn; EE(ei).YEK(:)];
end

% Sort points
[gi,ii] = sort(gi);
bi = bi(ii);
mi = mi(ii);
xn = xn(ii);
yn = yn(ii);

% Remove repeated nodes
[gi,ii,~] = unique(gi);
bi = logical(bi(ii));
mi = logical(mi(ii));
xn = xn(ii);
yn = yn(ii);

% Extract boundary/mortar nodes coordinates and indexes
%#ok<*FNDSB>
gb = gi(bi);
xb = xn(gb);
yb = yn(gb);

gm = gi(mi);
xm = xn(gm);
ym = yn(gm);

if verbose
  gcf
  plot(xb, yb, 'or','linewidth',2)
  plot(xm, ym, 'sg','linewidth',2)
  title('Mesh, nodes global indexes and boundary/mortar nodes');
  drawnow
end









%%%%%% Old version with E matrix and EE struct outputs
% Returns the mesh definition in two different formats. 
% E is a matrix containing, for each node, the following information:
% E = [elementIndex, nodes x coordinates, nodes y coordinates, ...
%      node GLL weights, node associated Jacobian, ...
%      node local index, node global index, node boundary flag, ...
%      node mortar flag];



% % Build element matrix from EE (E and EE should be unified at the end)
% E = [];
% for ei = 1 : length(EE)
%   E = [E; ...
%     ei*ones((G.N+1)^2,1), EE(ei).XEK(:), EE(ei).YEK(:), ...
%     EE(ei).AK(:), EE(ei).J(:), ...
%     EE(ei).li(:), EE(ei).gi(:), logical(EE(ei).bi(:)), logical(EE(ei).mi(:))];
%   
%   [~,I] = sort(E(:,7));
%   E     = E(I,:);
% end

