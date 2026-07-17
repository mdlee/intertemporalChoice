function path = resolveJagsModelFile(modelsDir, modelFileName)
%RESOLVEJAGSMODELFILE  Find a *_jags.txt in modelsDir or modelsDir/old.
%
%   PATH = RESOLVEJAGSMODELFILE(MODELSDIR, 'foo_jags.txt')
%   Prefers MODELSDIR, then MODELSDIR/old (archived non-entrop baselines).

if nargin < 2 || isempty(modelFileName)
  error('resolveJagsModelFile:badInput', 'Expected model file name.');
end
candidates = { ...
  fullfile(modelsDir, modelFileName), ...
  fullfile(modelsDir, 'old', modelFileName)};
for k = 1:numel(candidates)
  if isfile(candidates{k})
    path = candidates{k};
    return;
  end
end
error('resolveJagsModelFile:notFound', ...
  'JAGS model not found in %s or %s/old: %s', modelsDir, modelsDir, modelFileName);
end
