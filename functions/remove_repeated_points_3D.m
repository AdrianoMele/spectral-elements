function [x,y,u,idx] = remove_repeated_points_3D(xx,yy,uu)

idx = zeros(length(xx),1);
for i = 1 : length(xx)
%   temp = [xx(i) yy(i) uu(i)] == [xx(i+1:end) yy(i+1:end) uu(i+1:end)]; % uu can probably be removed
%   temp = temp(:,1)&temp(:,2)&temp(:,3);
  temp = [xx(i) yy(i)] == [xx(i+1:end) yy(i+1:end)]; % uu can probably be removed
  temp = temp(:,1)&temp(:,2);
  idx(i+find(temp)) = 1;
end
idx = ~idx;
x = xx(idx);
y = yy(idx);
u = uu(idx);