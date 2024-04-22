function [Xr,Yr,Ur] = refine_solution_SEM_2D(G,EE,u,gi,Nx,Ny,plotflag)
% [xr,yr,ur] = refine_solution_SEM_2D(EE,u,gi,Nx,Ny,[plotflag])
% Re-interpolate solution over a finer grid.
% EE is the mesh structure
% u is the solution associated to the GLL nodes
% gi are the global indexes of the nodes
% Nx,Ny are the number of points along the x/y direction per element
% plotflag: if =1, plots the solution
%
% Xr,Yr,Ur are three matrices containing the solution points.


% ue, xe, ye are matrices containing the coordinates of the interpolation
% points and the solution values for each element (each column is one
% element)


if nargin<7, plotflag = 0; end

% Global indexes may not be sorted
[~,ii] = sort(gi);
u = u(ii);

% Define interpolation points on reference element 
x = linspace(-1,1,Nx);
y = linspace(-1,1,Ny);
[Psi,XR,YR] = GLL_poly_2D(G.N,x,y);

% Initialize matrices
xe = zeros(Nx*Ny,length(EE));
ye = xe;
ue = xe;

% Cycle over elements
for ei = 1 : length(EE)
  UE = zeros(Ny,Nx);
  [XE,YE] = affine_transform(EE(ei).xve,EE(ei).yve,XR,YR);
  for i = 1 : G.N+1
    for j = 1 : G.N+1
      UE = UE + squeeze(Psi(i,j,:,:))*u(EE(ei).gi(i,j));
    end
  end
  xe(:,ei) = XE(:);
  ye(:,ei) = YE(:);
  ue(:,ei) = UE(:);
  
  % Plot solution
  if plotflag
    hold on
    mesh(XE,YE,UE);
  end
end

% Not efficient, but mesh-independent
% [xr,yr,ur] = remove_repeated_points_3D(xe(:),ye(:),ue(:)); % old
Nx = Nx*(G.MC-1)+1;
Ny = Ny*(G.MR-1)+1;
[Xr,Yr] = meshgrid(linspace(G.xmin,G.xmax,Nx),linspace(G.ymin,G.ymax,Ny));
warning off
Ur = griddata(xe,ye,ue,Xr,Yr);
warning on
if plotflag
  plot3(Xr,Yr,Ur,'.b','markersize',10)
end
