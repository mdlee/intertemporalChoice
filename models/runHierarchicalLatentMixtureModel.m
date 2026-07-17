function runHierarchicalLatentMixtureModel(modelName, initGenerator, varargin)
%RUNHIERARCHICALLATENTMIXTUREMODEL  Fit hierarchical latent mixture (all participants).
%
%   MCMC defaults: 12 chains (all required to converge),
%   burn-in 2e3 (fixed), 5e3 samples (fixed), thin 1 (×2 per failure).
%   preLoad (default true): use storage/{modelName}_*.mat only if max z R-hat
%   <= rhatCritical and enough chains pass; otherwise refit and overwrite.

p = inputParser;
addParameter(p, 'preLoad', true, @islogical);
addParameter(p, 'resetThin', false, @islogical);
addParameter(p, 'rhatCritical', 1.4, @isnumeric);
addParameter(p, 'keepChainsMin', 8, @isnumeric);
addParameter(p, 'nChains', 12, @isnumeric);
addParameter(p, 'nBurnin', 2e3, @isnumeric);
addParameter(p, 'nSamples', 5e3, @isnumeric);
addParameter(p, 'nThin', 1, @isnumeric);
addParameter(p, 'nModels', 11, @isnumeric);
addParameter(p, 'doParallel', true, @islogical);
addParameter(p, 'dataName', 'intertemporalChoice', @ischar);
addParameter(p, 'dataDir', '', @ischar);
addParameter(p, 'saveFigures', true, @islogical);
parse(p, varargin{:});

engine = 'jags';
modelsDir = fileparts(mfilename('fullpath'));
generalDir = fullfile(modelsDir, '..', 'general');
figuresDir = fullfile(modelsDir, 'figures');
storageDir = fullfile(modelsDir, 'storage');

addpath(generalDir);
cleanupObj = onCleanup(@() rmpath(generalDir));

[data, d] = prepareIntertemporalChoiceData(p.Results.dataName, p.Results.dataDir);
data.nModels = p.Results.nModels;

fileName = sprintf('%s_%s_%s.mat', modelName, p.Results.dataName, engine);
storagePath = fullfile(storageDir, fileName);

if p.Results.preLoad && isfile(storagePath)
  fprintf('Checking stored chains: %s\n', storagePath);
  S = load(storagePath, 'chains', 'stats', 'diagnostics', 'info');
  chains = S.chains;
  if isfield(S, 'stats'), stats = S.stats; end
  if isfield(S, 'diagnostics'), diagnostics = S.diagnostics; end
  if isfield(S, 'info'), info = S.info; end

  [zRhat, keepChains, rHatEach] = maxRhatOverParticipants( ...
    chains, p.Results.keepChainsMin, p.Results.rhatCritical);
  fprintf('z R-hat: max=%.3f, min=%.3f, %d/%d at 1 (constant z); keep %d / %d chains\n', ...
    zRhat, min(rHatEach), sum(rHatEach <= 1 + 1e-9), numel(rHatEach), ...
    numel(keepChains), p.Results.nChains);

  if zRhat <= p.Results.rhatCritical && numel(keepChains) >= p.Results.keepChainsMin
    fprintf('Stored chains meet R-hat <= %.3g — using saved fit\n', p.Results.rhatCritical);
    chains = subsetChainsLatentMixture(chains, keepChains);
    showFinalLatentMixtureFigures(chains, data, d, modelName, figuresDir, p.Results.saveFigures);
    return;
  end
  fprintf('Stored chains do not meet R-hat <= %.3g — refitting\n', p.Results.rhatCritical);
end

if ~isfolder(storageDir)
  mkdir(storageDir);
end
if p.Results.resetThin
  clearMcmcState(storageDir, modelName);
end
[nThin, nBurnin] = initialMcmcFromState(storageDir, modelName, p.Results.nThin, p.Results.nBurnin);
converged = false;

while ~converged
  tic;
  modelPath = resolveJagsModelFile(modelsDir, sprintf('%s_%s.txt', modelName, engine));
  [stats, chains, diagnostics, info] = callbayes(engine, ...
    'model', modelPath, ...
    'data', data, ...
    'outputname', 'samples', ...
    'init', initGenerator, ...
    'datafilename', modelName, ...
    'initfilename', modelName, ...
    'scriptfilename', modelName, ...
    'logfilename', fullfile('tmp', modelName), ...
    'nchains', p.Results.nChains, ...
    'nburnin', nBurnin, ...
    'nsamples', p.Results.nSamples, ...
    'monitorparams', {'z'}, ...
    'thin', nThin, ...
    'workingdir', fullfile('tmp', modelName), ...
    'verbosity', 0, ...
    'saveoutput', true, ...
    'allowunderscores', 1, ...
    'parallel', p.Results.doParallel);
  fprintf('%s took %.1f s (thin=%d, burn-in=%g)\n', upper(engine), toc, nThin, nBurnin);

  [zRhat, keepChains, rHatEach] = maxRhatOverParticipants( ...
    chains, p.Results.keepChainsMin, p.Results.rhatCritical);
  fprintf('z R-hat: max=%.3f, min=%.3f, %d/%d at 1 (constant z); keep %d / %d chains\n', ...
    zRhat, min(rHatEach), sum(rHatEach <= 1 + 1e-9), numel(rHatEach), ...
    numel(keepChains), p.Results.nChains);
  plotLatentMixtureZPosterior(chains, ...
    'nModels', data.nModels, ...
    'sgtitle', sprintf('%s | thin=%d | max R-hat=%.3f', modelName, nThin, zRhat), ...
    'savePath', hierarchicalFigurePath(figuresDir, modelName, sprintf('thin%d', nThin), p.Results.saveFigures));

  if zRhat <= p.Results.rhatCritical && numel(keepChains) >= p.Results.keepChainsMin
    chains = subsetChainsLatentMixture(chains, keepChains);
    converged = true;
    save(storagePath, 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
    fprintf('Saved %s\n', storagePath);
    clearMcmcState(storageDir, modelName);
  else
    saveUnsuccessfulMcmcState(storageDir, modelName, nThin);
    nThinNext = nextThin(nThin);
    fprintf('Not converged — next: thin=%d (was %d), burn-in=%g (unchanged)\n', ...
      nThinNext, nThin, nBurnin);
    nThin = nThinNext;
  end

  grtable(chains, p.Results.rhatCritical);
  codatable(chains);
end

showFinalLatentMixtureFigures(chains, data, d, modelName, figuresDir, p.Results.saveFigures);

end

function showFinalLatentMixtureFigures(chains, data, d, modelName, figuresDir, saveFigures)
zMean = get_matrix_from_coda(chains, 'z', @mean);
zMode = get_matrix_from_coda(chains, 'z', @mode);
disp(table((1:d.nParticipants)', zMean(:), zMode(:), 'VariableNames', {'p', 'zMean', 'zMode'}));
plotLatentMixtureZPosterior(chains, ...
  'nModels', data.nModels, ...
  'sgtitle', sprintf('%s | final', modelName), ...
  'savePath', hierarchicalFigurePath(figuresDir, modelName, 'final', saveFigures));
end

function chainsOut = subsetChainsLatentMixture(chainsIn, keepChains)
fn = fieldnames(chainsIn);
chainsOut = struct();
for k = 1:numel(fn)
  x = chainsIn.(fn{k});
  if isnumeric(x) && ndims(x) >= 2 && size(x, 2) >= max(keepChains)
    if ndims(x) == 2
      chainsOut.(fn{k}) = x(:, keepChains);
    else
      idx = repmat({':'}, 1, ndims(x));
      idx{2} = keepChains;
      chainsOut.(fn{k}) = x(idx{:});
    end
  else
    chainsOut.(fn{k}) = x;
  end
end
end
