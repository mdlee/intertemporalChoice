function [converged, attemptInfo] = runHierarchicalLatentMixtureEntropStan(modelName, varargin)
%RUNHIERARCHICALLATENTMIXTUREENTROPSTAN  Sequential Trinity Stan for the 11-component mixture.
%
%   Manuscript model: 8 hierarchical cognitive components + Guess/LL/SS.
%   Discrete z is marginalized; generated quantities draw z ~ categorical(softmax).
%   Gate (best subset of at least keepChainsMin of nChains):
%     max participant z R-hat <= rhatCritical on the trailing nSamplesMin,
%     and stacked n >= nSamplesMin.
%   Robustness jobs share the baseline .stan; pass 'muPrecScale' 0.5 or 2.

p = inputParser;
addParameter(p, 'preLoad', true, @islogical);
addParameter(p, 'rhatCritical', 1.10, @isnumeric);
addParameter(p, 'nSamplesMin', 1e4, @isnumeric);
addParameter(p, 'keepChainsMin', 4, @isnumeric);
addParameter(p, 'nChains', 6, @isnumeric);
addParameter(p, 'nBurnin', 1e3, @isnumeric);
addParameter(p, 'nBurninContinue', 2e2, @isnumeric);
addParameter(p, 'nSamples', 1e3, @isnumeric);
addParameter(p, 'nModels', 11, @isnumeric);
addParameter(p, 'doParallel', true, @islogical);
addParameter(p, 'dataName', 'intertemporalChoice', @ischar);
addParameter(p, 'dataDir', '', @ischar);
addParameter(p, 'saveFigures', true, @islogical);
addParameter(p, 'resetThin', false, @islogical);
addParameter(p, 'coldStart', false, @islogical);
addParameter(p, 'muPrecScale', 1, @isnumeric);
addParameter(p, 'saveOnConverge', true, @islogical);
addParameter(p, 'singleAttempt', false, @islogical);
parse(p, varargin{:});

modelsDir = fileparts(mfilename('fullpath'));
generalDir = fullfile(modelsDir, '..', 'general');
figuresDir = fullfile(modelsDir, 'stan', 'figures');
storageDir = fullfile(modelsDir, 'stan', 'storage');
logDir = fullfile(modelsDir, 'stan', 'logs');
if ~isfolder(logDir), mkdir(logDir); end
if ~isfolder(figuresDir), mkdir(figuresDir); end
if ~isfolder(storageDir), mkdir(storageDir); end
checkpointFile = fullfile(logDir, sprintf('%s.checkpoint.mat', modelName));
statusFile = fullfile(logDir, sprintf('%s.status', modelName));

addpath(generalDir);
cleanupObj = onCleanup(@() rmpath(generalDir)); %#ok<NASGU>

[data, d] = prepareIntertemporalChoiceData(p.Results.dataName, p.Results.dataDir);
data.nModels = p.Results.nModels;
data.muPrecScale = p.Results.muPrecScale;

nThin = 1;
nSamplesBatch = p.Results.nSamples;
nBurninFirst = p.Results.nBurnin;
nBurninContinue = p.Results.nBurninContinue;
nSamplesMin = p.Results.nSamplesMin;
checkpointKeep = max(nSamplesMin, 2 * nSamplesBatch);
nChains = p.Results.nChains;
nP = double(d.nParticipants);
nModels = p.Results.nModels;

rawInit = makeLatentMixtureEntropInit(storageDir, p.Results.dataName, 'stan', ...
  nP, 'nModels', nModels, 'jitterFrac', 0.05, 'seed', 1);
initGenerator = @() dropNonParameters(rawInit());
monitorParams = mixtureEntropStanMonitorParams();

storagePath = storageMatPath(storageDir, modelName, p.Results.dataName, 'stan', '');

attemptInfo = struct('rMax', nan, 'keepChains', [], 'nCollected', 0, ...
  'batch', 0, 'converged', false);

if p.Results.preLoad && isfile(storagePath)
  fprintf('Checking stored chains: %s\n', storagePath);
  S = load(storagePath, 'chains');
  chains = S.chains;
  nStored = codaNSamples(chains);
  gateChains = tailCoda(chains, min(nStored, nSamplesMin));
  [ok, keepChains, rMax] = findConvergedZChainSubset( ...
    gateChains, p.Results.keepChainsMin, p.Results.rhatCritical);
  attemptInfo.rMax = rMax;
  attemptInfo.keepChains = keepChains;
  attemptInfo.nCollected = nStored;
  fprintf('Stored z R-hat=%.3f keep %d/%d n=%d\n', rMax, numel(keepChains), nChains, nStored);
  if ok && nStored >= nSamplesMin
    fprintf('Stored mixture meets gate — using saved fit\n');
    plotLatentMixtureZPosterior(chains, 'nModels', nModels, ...
      'sgtitle', sprintf('%s | stored', modelName), ...
      'savePath', hierarchicalFigurePath(figuresDir, modelName, 'final', p.Results.saveFigures));
    attemptInfo.converged = true;
    converged = true;
    return;
  end
end

if p.Results.resetThin || p.Results.coldStart
  if isfile(checkpointFile)
    delete(checkpointFile);
  end
end

[acc, lastInitCell, initFn, nBurnin, bestRhat, nCollected, batch] = ...
  resumeMixtureStanCheckpoint(checkpointFile, modelName, initGenerator, ...
    nBurninFirst, nBurninContinue);

converged = false;
stats = struct();
diagnostics = struct();
info = struct();

while true
  batch = batch + 1;
  fprintf('\n--- %s | batch=%d | burnin=%g | collect=%g | stacked=%g ---\n', ...
    modelName, batch, nBurnin, nSamplesBatch, nCollected);
  fprintf('Wall time %s\n', datestr(now, 31));
  writeMcmcStatus(statusFile, sprintf( ...
    '%s | batch=%d | stacked=%g | collect +%g | best z R-hat=%.3f | running', ...
    modelName, batch, nCollected, nSamplesBatch, bestRhat));

  [stats, chains, diagnostics, info] = fitMixtureStanOnce( ...
    modelsDir, modelName, monitorParams, data, initFn, ...
    nThin, nBurnin, nSamplesBatch, nChains, p.Results.doParallel);

  acc = appendCoda(acc, chains);
  acc = trimCodaTail(acc, checkpointKeep);
  if batch == 1 && nCollected == 0
    nCollected = codaNSamples(acc);
  else
    nCollected = nCollected + nSamplesBatch;
  end
  lastInitCell = codaLastInits(chains, nP, nChains, nModels);
  lastInitCell = cellfun(@dropNonParameters, lastInitCell, 'UniformOutput', false);
  initFn = makeSerialInitFn(lastInitCell);
  nBurnin = nBurninContinue;

  gateAcc = tailCoda(acc, min(codaNSamples(acc), nSamplesMin));
  nGate = codaNSamples(gateAcc);
  [ok, keepChains, rMax, rHatEach] = findConvergedZChainSubset( ...
    gateAcc, p.Results.keepChainsMin, p.Results.rhatCritical);
  fprintf(['stacked n=%d | trail n=%d | z R-hat max=%.3f | keep %s (%d/%d) | ' ...
    '%d/%d participants at 1\n'], ...
    nCollected, nGate, rMax, mat2str(keepChains), numel(keepChains), nChains, ...
    sum(rHatEach <= 1 + 1e-9), numel(rHatEach));

  try
    plotLatentMixtureZPosterior(gateAcc, 'nModels', nModels, ...
      'sgtitle', sprintf('%s | trail=%d | z R-hat=%.3f', modelName, nGate, rMax), ...
      'savePath', hierarchicalFigurePath(figuresDir, modelName, ...
        sprintf('trail%d', nGate), p.Results.saveFigures));
  catch ME
    fprintf('z figure skipped: %s\n', ME.message);
  end

  if rMax < bestRhat
    bestRhat = rMax;
    if p.Results.saveOnConverge
      saveMixtureBest(storagePath, gateAcc, keepChains, p.Results.keepChainsMin, ...
        stats, diagnostics, info);
      fprintf('Updated best-so-far (z R-hat=%.3f, n=%d) → %s\n', bestRhat, nGate, storagePath);
    end
  end

  saveSeqCheckpoint(checkpointFile, modelName, acc, lastInitCell, ...
    bestRhat, nCollected, batch);

  converged = ok && nCollected >= nSamplesMin;
  attemptInfo.rMax = rMax;
  attemptInfo.keepChains = keepChains;
  attemptInfo.nCollected = nCollected;
  attemptInfo.batch = batch;
  attemptInfo.converged = converged;
  writeMcmcStatus(statusFile, sprintf( ...
    '%s | batch=%d | stacked=%d | trail=%d | z R-hat=%.3f | best=%.3f | keep=%d/%d | converged=%d', ...
    modelName, batch, nCollected, nGate, rMax, bestRhat, numel(keepChains), ...
    nChains, converged));

  if converged || p.Results.singleAttempt
    break;
  end
  fprintf('Appending another %g samples\n', nSamplesBatch);
end

if converged
  chains = subsetChainsByKeep(gateAcc, keepChains);
  if p.Results.saveOnConverge
    save(storagePath, 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
    fprintf('Saved %s (%d trailing samples)\n', storagePath, codaNSamples(chains));
    if isfile(checkpointFile)
      delete(checkpointFile);
    end
    plotLatentMixtureZPosterior(chains, 'nModels', nModels, ...
      'sgtitle', sprintf('%s | final', modelName), ...
      'savePath', hierarchicalFigurePath(figuresDir, modelName, 'final', p.Results.saveFigures));
  end
end
end

function names = mixtureEntropStanMonitorParams()
names = { ...
  'z', 'w', 'muW', 'tauW', ...
  'kappaEX', 'muKappaEX', 'tauKappaEX', ...
  'kappaHC', 'muKappaHC', 'tauKappaHC', ...
  'kappaHY', 'tauHY', 'muKappaHY', 'tauKappaHY', 'muTauHY', 'tauTauHY', ...
  'deltaPD', 'muDeltaPD', 'tauDeltaPD', ...
  'deltaDD', 'omegaDD', 'muDeltaDD', 'tauDeltaDD', 'muOmegaDD', 'tauOmegaDD', ...
  'gammaTR', 'tauTR', 'kappaTR', 'varthetaTR', ...
  'muGammaTR', 'tauGammaTR', 'muTauTR', 'tauTauTR', 'muKappaTR', 'tauKappaTR', ...
  'muVarthetaTR', 'tauVarthetaTR', ...
  'gammaUT', 'tauUT', 'kappaUT', 'varthetaUT', 'etaUT', ...
  'muGammaUT', 'tauGammaUT', 'muTauUT', 'tauTauUT', 'muKappaUT', 'tauKappaUT', ...
  'muVarthetaUT', 'tauVarthetaUT', 'muEtaUT', 'tauEtaUT', ...
  'beta0IT', 'betaRAIT', 'betaRRIT', 'betaTAIT', 'betaTRIT', ...
  'muBeta0IT', 'tauBeta0IT', 'muBetaRAIT', 'tauBetaRAIT', 'muBetaRRIT', 'tauBetaRRIT', ...
  'muBetaTAIT', 'tauBetaTAIT', 'muBetaTRIT', 'tauBetaTRIT', ...
  'alphaLL', 'alphaSS'};
end

function [stats, chains, diagnostics, info] = fitMixtureStanOnce( ...
    modelsDir, modelName, monitorParams, data, initFn, ...
    nThin, nBurnin, nSamples, nChains, doParallel)
tic;
[modelPath, workDir] = stageStanEntropModel(modelsDir, modelName);
[stats, chains, diagnostics, info] = callbayes('stan', ...
  'model', modelPath, ...
  'data', data, ...
  'outputname', 'samples', ...
  'init', initFn, ...
  'datafilename', modelName, ...
  'initfilename', modelName, ...
  'scriptfilename', modelName, ...
  'logfilename', fullfile(workDir, 'log'), ...
  'nchains', nChains, ...
  'nburnin', nBurnin, ...
  'nsamples', nSamples, ...
  'monitorparams', monitorParams, ...
  'thin', nThin, ...
  'workingdir', workDir, ...
  'verbosity', 0, ...
  'saveoutput', true, ...
  'allowunderscores', 1, ...
  'parallel', doParallel);
fprintf('STAN took %.1f s (thin=%d, burn-in=%g, collect=%g)\n', ...
  toc, nThin, nBurnin, nSamples);
end

function s = dropNonParameters(s)
drop = {'z', 'lp__', 'accept_stat__', 'stepsize__', 'treedepth__', ...
  'n_leapfrog__', 'divergent__', 'energy__', 'deviance'};
for k = 1:numel(drop)
  if isfield(s, drop{k})
    s = rmfield(s, drop{k});
  end
end
end

function acc = appendCoda(acc, new)
if isempty(acc)
  acc = new;
  return;
end
fn = fieldnames(new);
for k = 1:numel(fn)
  name = fn{k};
  b = new.(name);
  if isfield(acc, name)
    a = acc.(name);
    if isnumeric(a) && isnumeric(b) && size(a, 2) == size(b, 2)
      acc.(name) = [a; b];
    else
      acc.(name) = b;
    end
  else
    acc.(name) = b;
  end
end
end

function acc = trimCodaTail(acc, keep)
n = codaNSamples(acc);
if n <= keep
  return;
end
fn = fieldnames(acc);
for k = 1:numel(fn)
  x = acc.(fn{k});
  if isnumeric(x) && size(x, 1) == n
    acc.(fn{k}) = x(end - keep + 1:end, :);
  end
end
end

function t = tailCoda(acc, nKeep)
n = codaNSamples(acc);
if n <= nKeep
  t = acc;
  return;
end
t = acc;
fn = fieldnames(acc);
for k = 1:numel(fn)
  x = acc.(fn{k});
  if isnumeric(x) && size(x, 1) == n
    t.(fn{k}) = x(end - nKeep + 1:end, :);
  end
end
end

function n = codaNSamples(chains)
if isempty(chains)
  n = 0;
  return;
end
fn = fieldnames(chains);
n = 0;
for k = 1:numel(fn)
  x = chains.(fn{k});
  if isnumeric(x) && ndims(x) == 2
    n = size(x, 1);
    return;
  end
end
end

function initCell = codaLastInits(chains, nP, nChains, nModels)
fn = fieldnames(chains);
initCell = cell(1, nChains);
scalarNames = {};
indexed = struct();
skip = {'lp__', 'accept_stat__', 'stepsize__', 'treedepth__', ...
  'n_leapfrog__', 'divergent__', 'energy__', 'deviance'};
for k = 1:numel(fn)
  name = fn{k};
  if any(strcmp(name, skip))
    continue;
  end
  tok = regexp(name, '^([A-Za-z]\w*)_(\d+)$', 'tokens', 'once');
  if isempty(tok)
    scalarNames{end + 1} = name; %#ok<AGROW>
  else
    base = tok{1};
    idx = str2double(tok{2});
    if ~isfield(indexed, base)
      indexed.(base) = [];
    end
    indexed.(base)(end + 1) = idx; %#ok<AGROW>
  end
end
bases = fieldnames(indexed);
for c = 1:nChains
  s = struct();
  for b = 1:numel(bases)
    base = bases{b};
    if strcmp(base, 'z')
      continue;
    end
    ids = sort(indexed.(base));
    v = zeros(max(ids), 1);
    for ii = 1:numel(ids)
      fld = sprintf('%s_%d', base, ids(ii));
      v(ids(ii)) = chains.(fld)(end, c);
    end
    if numel(v) >= nP
      v = v(1:nP);
    end
    s.(base) = v(:);
  end
  for b = 1:numel(scalarNames)
    nm = scalarNames{b};
    s.(nm) = chains.(nm)(end, c);
  end
  initCell{c} = s;
end
end

function fh = makeSerialInitFn(initCell)
idx = 0;
fh = @nextInit;
  function s = nextInit
    idx = idx + 1;
    if idx > numel(initCell)
      error('makeSerialInitFn:exhausted', 'Init list exhausted at %d.', idx);
    end
    s = initCell{idx};
  end
end

function saveMixtureBest(path, acc, keepChains, keepChainsMin, stats, diagnostics, info)
if numel(keepChains) >= keepChainsMin
  chains = subsetChainsByKeep(acc, keepChains);
else
  chains = acc;
end
save(path, 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
end

function chainsOut = subsetChainsByKeep(chainsIn, keepChains)
fn = fieldnames(chainsIn);
chainsOut = struct();
for k = 1:numel(fn)
  x = chainsIn.(fn{k});
  if isnumeric(x) && ndims(x) >= 2 && size(x, 2) >= max(keepChains)
    chainsOut.(fn{k}) = x(:, keepChains);
  else
    chainsOut.(fn{k}) = x;
  end
end
end

function saveSeqCheckpoint(path, modelName, acc, lastInitCell, bestRhat, nCollected, batch)
save(path, 'modelName', 'acc', 'lastInitCell', 'bestRhat', 'nCollected', 'batch', '-v7.3');
end

function [acc, lastInitCell, initFn, nBurnin, bestRhat, nCollected, batch] = ...
    resumeMixtureStanCheckpoint(checkpointFile, modelName, initGenerator, ...
    nBurninFirst, nBurninContinue)
acc = [];
lastInitCell = {};
initFn = initGenerator;
nBurnin = nBurninFirst;
bestRhat = inf;
nCollected = 0;
batch = 0;
if ~isfile(checkpointFile)
  return;
end
ck = load(checkpointFile);
if isfield(ck, 'modelName') && strcmp(ck.modelName, modelName) ...
    && isfield(ck, 'lastInitCell') && ~isempty(ck.lastInitCell)
  acc = ck.acc;
  lastInitCell = ck.lastInitCell;
  initFn = makeSerialInitFn(lastInitCell);
  nBurnin = nBurninContinue;
  bestRhat = ck.bestRhat;
  nCollected = ck.nCollected;
  batch = ck.batch;
  fprintf('Resuming mixture Stan checkpoint: stacked=%d batch=%d best z R-hat=%.3f\n', ...
    nCollected, batch, bestRhat);
end
end
