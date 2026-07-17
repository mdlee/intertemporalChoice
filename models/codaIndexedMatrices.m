function matrices = codaIndexedMatrices(coda, paramName)
%CODAINDEXEDMATRICES  Trinity CODA fields param_1, param_2, ... as sample matrices.
%
%   MATRICES = CODAINDEXEDMATRICES(CODA) with default paramName 'z' returns a
%   cell array of nSamples-by-nChains matrices, one per index. Also accepts a
%   single matrix CODA.(paramName) for scalar parameters.

if nargin < 2 || isempty(paramName)
  paramName = 'z';
end

if isnumeric(coda)
  matrices = {coda};
  return;
end

if ~isstruct(coda)
  error('codaIndexedMatrices:badInput', 'Expected struct or numeric matrix.');
end

if isfield(coda, paramName) && isstruct(coda.(paramName))
  coda = coda.(paramName);
end

if isfield(coda, paramName) && isnumeric(coda.(paramName))
  matrices = {coda.(paramName)};
  return;
end

[selection, nSel] = trinity.select_fields(coda, ['^' paramName '_']);
if nSel == 0
  error('codaIndexedMatrices:noParam', ...
    'No fields matching %s_1, %s_2, ... in CODA.', paramName, paramName);
end

idx = zeros(nSel, 1);
for k = 1:nSel
  tail = selection{k}(numel(paramName) + 2:end);
  idx(k) = str2double(tail);
  if isnan(idx(k))
    error('codaIndexedMatrices:badField', 'Cannot parse index from %s.', selection{k});
  end
end

[~, ord] = sort(idx);
selection = selection(ord);
matrices = cell(nSel, 1);
for k = 1:nSel
  x = coda.(selection{k});
  if ~isnumeric(x) || ndims(x) ~= 2
    error('codaIndexedMatrices:badField', '%s is not nSamples-by-nChains.', selection{k});
  end
  matrices{k} = x;
end

end
