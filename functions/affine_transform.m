function [xe,ye,J,M,Jinv] = affine_transform(xve,yve,xr,yr)
% [xe,ye] = affine_transform(xve,xye,xr,yr)
% Returns the position of the points (xr,yr) defined in the reference
% element [-1,1]x[-1,1] when transported on the element of vertexes 
% (xve,yve).
%
% [xe,ye,J] = affine_transform(xve,xye,xr,yr) 
% also returns the Jacobian of the transformation [xe(xr,yr), ye(xr,yr)] at
% each point. 
%
% [xe,ye,J,M] = affine_transform(xve,xye,xr,yr) 
% also returns measure of the element M.
%
% [xe,ye,J,M,Jinv] = affine_transform(xve,xye,xr,yr) 
% also the full inverse Jacobian matrix in each point. The dimension of JJ
% is [2]x[2]x[Nx]x[Ny].
%
% Th vertexes of the reference elements are:
%    xvr = [-1 +1 +1 -1];
%    yvr = [-1 -1 +1 +1];
%
% The mapping is obtained as:
%    xe = a1*xr + b1*yr + c1*xr*yr + d1;
%    xe = a2*xr + b2*yr + c2*xr*yr + d2;
%
% The transform matrix is obtained as the inverse of:
%    A = [-1 -1  +1  1   0  0   0  0
%         +1 -1  -1  1   0  0   0  0
%         +1 +1  +1  1   0  0   0  0
%         -1 +1  -1  1   0  0   0  0
%          0  0   0  0  -1 -1  +1  1
%          0  0   0  0  +1 -1  -1  1
%          0  0   0  0  +1 +1  +1  1
%          0  0   0  0  -1 +1  -1  1];
% i.e. A\[xe ye]' = [a1 b1 c1 d1 a2 b2 c2 d2]';
%
% The elements of the Jacobian matrix are:
%     J(1,1) = dxe/dxr = a1+c1*yr;
%     J(1,2) = dxe/dyr = b1+c1*xr;
%     J(2,1) = dye/dxr = a2+c2*yr;
%     J(2,2) = dye/dyr = b1+c2*xr;
% i.e. J*[dxr; dyr] = [dxe; dye]. The inverse Jacobian is obtained as 
% inv(J).


if iscolumn(xve), xve = xve'; end
if iscolumn(yve), yve = yve'; end

% % Sort vertexes anticlockwise
% tve = atan2(yve-mean(xve),xve-mean(yve));
% [~, idx] = sort(tve);
% xve = xve(idx); yve = yve(idx);

% Compute coefficients of the transformation
b = [xve(1:4), yve(1:4)]';
   
A = [-0.25    0.25    0.25   -0.25       0       0       0       0
     -0.25   -0.25    0.25    0.25       0       0       0       0
      0.25   -0.25    0.25   -0.25       0       0       0       0
      0.25    0.25    0.25    0.25       0       0       0       0
         0       0       0       0   -0.25    0.25    0.25   -0.25
         0       0       0       0   -0.25   -0.25    0.25    0.25
         0       0       0       0    0.25   -0.25    0.25   -0.25
         0       0       0       0    0.25    0.25    0.25    0.25];
       
coef = A*b;

B = [coef([1 2 3]) coef([5 6 7])]'; % [a1 b1 c1; a2 b2 c2]'
c = coef([4 8]);                    % [d1 d2]

% Compute Jacobian in the points (xr,yr)
xe = xr*0; ye = xe; J = xr*0;
for i = 1 : size(xr,1)
  for j = 1 : size(xr,2)
    xy = B*[xr(i,j); yr(i,j); xr(i,j)*yr(i,j)] + c;
    xe(i,j) = xy(1); ye(i,j) = xy(2);
    
    J(i,j) = (B(1,1)+B(1,3)*xy(2))*(B(2,2)+B(2,3)*xy(1)) - ...
             (B(1,2)+B(1,3)*xy(1))*(B(2,1)+B(2,3)*xy(2));         
  end
end

% Compute polygon area
if nargout >= 4
  M = polyarea(xve,yve);
end

% Compute full inverse Jacobian matrix
Jinv = zeros(2,2,size(xr,1),size(xr,2));
if nargout >= 5
  for i = 1 : size(xr,1)
    for j = 1 : size(xr,2)     
      Jinv(:,:,i,j) = inv([(B(1,1)+B(1,3)*xy(2)) (B(1,2)+B(1,3)*xy(1))
                           (B(2,1)+B(2,3)*xy(2)) (B(2,2)+B(2,3)*xy(1))]);
    end
  end
end