function path = storageMatPath(storageDir, modelName, dataName, engine, storageTag)
%STORAGEMATPATH  Build storage .mat path with optional tag (does not overwrite).
%
%   Canonical (drawFiguresEntrop default):
%     {modelName}_{dataName}_{engine}.mat
%   Tagged archive / alternate run:
%     {modelName}_{dataName}_{engine}__{storageTag}.mat
%
%   STORAGETAG empty / omitted → canonical path.

if nargin < 5 || isempty(storageTag)
  storageTag = '';
end
storageTag = char(storageTag);
if isempty(storageTag)
  fileName = sprintf('%s_%s_%s.mat', modelName, dataName, engine);
else
  fileName = sprintf('%s_%s_%s__%s.mat', modelName, dataName, engine, storageTag);
end
path = fullfile(storageDir, fileName);
end
