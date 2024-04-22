function Deltapsi = laplacian(R,Z,Psi)
% (nabla2)psi = (d^2/dr^2 + d^2/dz^2) * psi

dr = R(1,2)-R(1,1);
dz = Z(1,1)-Z(2,1);
nr = size(R,2);
nz = size(Z,1);

% To be pre-multiplied
ddz = zeros(nz,nz);
for i = 2 : nz-1
    ddz(i,i-1) = -1;
    ddz(i,i+1) = +1;
end
ddz = ddz/2/dz;

% To be post-multiplied
ddr = zeros(nr,nr);
for j = 2 : nr-1
    ddr(j-1,j) = -1;
    ddr(j+1,j) = +1;
end
ddr = ddr/2/dr;

% Deltapsi = ddz*ddz*Psi + (((Psi*ddr)./R)*ddr).*R;
Deltapsi = ddz*ddz*Psi + Psi*ddr*ddr;

Deltapsi(1,:) = Deltapsi(3,:);
Deltapsi(:,1) = Deltapsi(:,3);
Deltapsi(2,:) = Deltapsi(3,:);
Deltapsi(:,2) = Deltapsi(:,3);
Deltapsi(end,:) = Deltapsi(end-2,:);
Deltapsi(:,end) = Deltapsi(:,end-2);
Deltapsi(end-1,:) = Deltapsi(end-2,:);
Deltapsi(:,end-1) = Deltapsi(:,end-2);









