function path = resolveJagsModelFile(modelsDir, modelFileName)
%RESOLVEJAGSMODELFILE  Find a *_jags.txt under models/.
%
%   PATH = RESOLVEJAGSMODELFILE(MODELSDIR, 'foo_jags.txt')
%   Search order:
%     modelsDir/jags/       — main models (original priors)
%     modelsDir/jagsRobust/ — half/double prior robustness variants
%     modelsDir/            — legacy flat layout
%     modelsDir/old/        — archived baselines

if nargin < 2 || isempty(modelFileName)
  error('resolveJagsModelFile:badInput', 'Expected model file name.');
end
candidates = { ...
  fullfile(modelsDir, 'jags', modelFileName), ...
  fullfile(modelsDir, 'jagsRobust', modelFileName), ...
  fullfile(modelsDir, modelFileName), ...
  fullfile(modelsDir, 'old', modelFileName)};
for k = 1:numel(candidates)
  if isfile(candidates{k})
    path = candidates{k};
    return;
  end
end
error('resolveJagsModelFile:notFound', ...
  'JAGS model not found for %s (searched jags/, jagsRobust/, ., old/ under %s)', ...
  modelFileName, modelsDir);
end
