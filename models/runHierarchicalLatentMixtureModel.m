function [converged, attemptInfo] = runHierarchicalLatentMixtureModel(modelName, initGenerator, varargin)
%RUNHIERARCHICALLATENTMIXTUREMODEL  Fit hierarchical latent mixture (all participants).
%
%   MCMC defaults: 12 chains (all required to converge),
%   burn-in 2e3 (fixed), 5e3 samples (fixed), thin 1 (×2 per failure).
%   preLoad (default true): use storage/{modelName}_*.mat only if max z R-hat
%   <= rhatCritical and enough chains pass; otherwise refit and overwrite.
%
%   Name-value extras:
%     storageTag     — if non-empty, save/load
%                      {model}_{data}_{engine}__{tag}.mat (canonical untouched)
%     singleAttempt  — one MCMC attempt at nThin; do not double thin on failure
%     saveOnConverge — save when gate passes (default true)
%     saveBestPath   — if non-empty, always write here when this attempt beats
%                      previousBestRhat (for drawFiguresEntrop canonical path)
%     previousBestRhat — scalar; used with saveBestPath (default Inf)

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
addParameter(p, 'storageTag', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'singleAttempt', false, @islogical);
addParameter(p, 'saveOnConverge', true, @islogical);
addParameter(p, 'saveBestPath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'previousBestRhat', inf, @isnumeric);
parse(p, varargin{:});

engine = 'jags';
modelsDir = fileparts(mfilename('fullpath'));
generalDir = fullfile(modelsDir, '..', 'general');
figuresDir = fullfile(modelsDir, 'figures');
storageDir = fullfile(modelsDir, 'storage');
storageTag = char(p.Results.storageTag);
saveBestPath = char(p.Results.saveBestPath);

addpath(generalDir);
cleanupObj = onCleanup(@() rmpath(generalDir));

[data, d] = prepareIntertemporalChoiceData(p.Results.dataName, p.Results.dataDir);
data.nModels = p.Results.nModels;

storagePath = storageMatPath(storageDir, modelName, p.Results.dataName, engine, storageTag);

attemptInfo = struct( ...
  'rMax', nan, ...
  'rHatEach', [], ...
  'nThin', p.Results.nThin, ...
  'keepChains', [], ...
  'chains', [], ...
  'bestRhat', p.Results.previousBestRhat, ...
  'savedBest', false, ...
  'converged', false);

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

  attemptInfo.rMax = zRhat;
  attemptInfo.rHatEach = rHatEach;
  attemptInfo.keepChains = keepChains;
  attemptInfo.chains = chains;
  attemptInfo.bestRhat = min(attemptInfo.bestRhat, zRhat);

  if zRhat <= p.Results.rhatCritical && numel(keepChains) >= p.Results.keepChainsMin
    fprintf('Stored chains meet R-hat <= %.3g — using saved fit\n', p.Results.rhatCritical);
    chains = subsetChainsLatentMixture(chains, keepChains);
    attemptInfo.chains = chains;
    attemptInfo.converged = true;
    showFinalLatentMixtureFigures(chains, data, d, modelName, figuresDir, p.Results.saveFigures);
    converged = true;
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
if p.Results.singleAttempt
  nThin = p.Results.nThin;
  nBurnin = p.Results.nBurnin;
else
  [nThin, nBurnin] = initialMcmcFromState(storageDir, modelName, p.Results.nThin, p.Results.nBurnin);
end
converged = false;
bestRhat = p.Results.previousBestRhat;

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

  attemptInfo.rMax = zRhat;
  attemptInfo.rHatEach = rHatEach;
  attemptInfo.nThin = nThin;
  attemptInfo.keepChains = keepChains;
  attemptInfo.chains = chains;

  if zRhat < bestRhat
    bestRhat = zRhat;
    attemptInfo.bestRhat = bestRhat;
    if ~isempty(saveBestPath)
      if numel(keepChains) >= p.Results.keepChainsMin
        chainsSave = subsetChainsLatentMixture(chains, keepChains);
      else
        chainsSave = chains;
      end
      toSave = struct('chains', chainsSave, 'stats', stats, ...
        'diagnostics', diagnostics, 'info', info);
      save(saveBestPath, '-struct', 'toSave', '-v7.3');
      attemptInfo.savedBest = true;
      fprintf('Updated best-so-far (z R-hat=%.3f) → %s\n', bestRhat, saveBestPath);
    end
  end

  if zRhat <= p.Results.rhatCritical && numel(keepChains) >= p.Results.keepChainsMin
    chains = subsetChainsLatentMixture(chains, keepChains);
    attemptInfo.chains = chains;
    converged = true;
    attemptInfo.converged = true;
    if p.Results.saveOnConverge
      save(storagePath, 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
      fprintf('Saved %s\n', storagePath);
      clearMcmcState(storageDir, modelName);
    end
  else
    saveUnsuccessfulMcmcState(storageDir, modelName, nThin);
    if p.Results.singleAttempt
      fprintf('Not converged at thin=%d (max z R-hat=%.3f)\n', nThin, zRhat);
      break;
    end
    nThinNext = nextThin(nThin);
    fprintf('Not converged — next: thin=%d (was %d), burn-in=%g (unchanged)\n', ...
      nThinNext, nThin, nBurnin);
    nThin = nThinNext;
  end

  grtable(chains, p.Results.rhatCritical);
  codatable(chains);

  if converged || p.Results.singleAttempt
    break;
  end
end

if converged || (~isempty(attemptInfo.chains) && isstruct(attemptInfo.chains))
  showFinalLatentMixtureFigures(attemptInfo.chains, data, d, modelName, figuresDir, p.Results.saveFigures);
end

end

function showFinalLatentMixtureFigures(chains, data, d, modelName, figuresDir, saveFigures)
if isempty(chains) || ~isstruct(chains)
  return;
end
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
