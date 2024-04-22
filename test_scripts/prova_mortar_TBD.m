lm = EE(ei).li(EE(ei).mi); %  local index of the mortar points
gm = EE(ei).gi(EE(ei).mi); % global index of the mortar points

% Extract information from EE
AAA = blkdiag(EE(:).AA);
FFF = vertcat(EE(:).FF);
gg  = vertcat(EE(:).GI);
mm  = vertcat(EE(:).MI); % ismortar
bb  = vertcat(EE(:).BI); % isboundary

% Remove boundary nodes
AAA = AAA(~bb,~bb);
FFF = FFF(~bb);
gg  = gg(~bb);
mm  = mm(~bb);

% Apply continuity at interface
[gn, i1, i2] = unique(gg);
S = sparse(diag(mm(i1)));
S = zeros(length(i1),length(i2));
for i = 1:length(i1)
  for j = 1:length(i2)
    S(i,j) = gg(i1(i))==gg(i2(j));
  end
end

ut = AAA\FFF;
Un(gin) = ut(i1);
Un(gb) = ub;