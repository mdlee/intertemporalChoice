function Z = loadMixtureZDraws(matPath, nParticipants)
%LOADMIXTUREZDRAWS  nDraws-by-nParticipants matrix of stored mixture z samples.

if ~isfile(matPath)
  error('loadMixtureZDraws:missing', 'Missing mixture mat: %s', matPath);
end
S = load(matPath, 'chains');
Z = [];
for i = 1:nParticipants
  fld = sprintf('z_%d', i);
  if ~isfield(S.chains, fld)
    error('loadMixtureZDraws:noZ', 'No %s in %s', fld, matPath);
  end
  col = S.chains.(fld)(:);
  if isempty(Z)
    Z = zeros(numel(col), nParticipants);
  elseif numel(col) ~= size(Z, 1)
    n = min(numel(col), size(Z, 1));
    Z = Z(1:n, :);
    col = col(1:n);
  end
  Z(:, i) = col;
end
end
