function [converged, attemptInfo] = runHierarchicalExecutionModel(modelName, monitorParams, initGenerator, varargin)
%RUNHIERARCHICALEXECUTIONMODEL  Fit one hierarchical execution model (all participants).
%
%   Sequential batches (same schedule as the latent mixture):
%     thin = 1, 6 chains,
%     first batch burn-in 1e3 then collect 1e3,
%     later batches: last-state re-init, burn-in 200, collect 1e3, append.
%   Continue until some subset of >= keepChainsMin chains has max
%   monitored-parameter R-hat <= rhatCritical (default 1.05) on the trailing
%   window of nSamplesMin (default 1e4) samples per chain. Sequential overnight
%   uses keepChainsMin = 4 (same subset idea as mixture robustness).
%   and total nCollected >= nSamplesMin.
%   Checkpoints retain only the last checkpointKeep samples (default 2*nSamplesMin).
%   preLoad (default true): skip refit if trailing-window gate passes on storage.
%   coldStart (default false): ignore canonical storage seeding (fresh chains).
%   seedStoragePath (''): optional .mat to seed checkpoint state (used even when
%     coldStart); overrides canonical storage for initialization only.
%   Checkpoint: models/logs/{modelName}.checkpoint.mat
%   rhatParticipants optionally restricts the convergence gate to selected
%   participants; all participants remain in the hierarchical model.
%
%   [CONVERGED, ATTEMPTINFO] = ... also returns rMax, rHatByParam, nThin, keepChains.

p = inputParser;
addParameter(p, 'preLoad', true, @islogical);
addParameter(p, 'singleAttempt', false, @islogical);
addParameter(p, 'figuresOnly', false, @islogical);
addParameter(p, 'plotFailedRuns', true, @islogical);
addParameter(p, 'rhatCritical', 1.05, @isnumeric);
addParameter(p, 'nSamplesMin', 1e4, @isnumeric);
addParameter(p, 'keepChainsMin', 6, @isnumeric);
addParameter(p, 'nChains', 6, @isnumeric);
addParameter(p, 'nBurnin', 1e3, @isnumeric);
addParameter(p, 'nBurninContinue', 2e2, @isnumeric);
addParameter(p, 'nSamples', 1e3, @isnumeric);
addParameter(p, 'nThin', 1, @isnumeric);
addParameter(p, 'doParallel', true, @islogical);
addParameter(p, 'dataName', 'intertemporalChoice', @ischar);
addParameter(p, 'dataDir', '', @ischar);
addParameter(p, 'saveFigures', true, @islogical);
addParameter(p, 'resetThin', false, @islogical);
addParameter(p, 'coldStart', false, @islogical);
addParameter(p, 'seedStoragePath', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'saveOnConverge', true, @islogical);
addParameter(p, 'checkpointKeep', [], @isnumeric);
addParameter(p, 'rhatParticipants', [], @isnumeric);
addParameter(p, 'engine', 'jags', @(s) ischar(s) || isstring(s));
addParameter(p, 'muPrecScale', 1, @isnumeric);
parse(p, varargin{:});

nThin = 1;
nSamplesBatch = p.Results.nSamples;
nBurninFirst = p.Results.nBurnin;
nBurninContinue = p.Results.nBurninContinue;
nSamplesMin = p.Results.nSamplesMin;
checkpointKeep = p.Results.checkpointKeep;
if isempty(checkpointKeep)
  checkpointKeep = max(nSamplesMin, 2 * nSamplesBatch);
end
if p.Results.nThin ~= 1
  fprintf('Ignoring nThin=%g — sequential batches always use thin=1\n', p.Results.nThin);
end

attemptInfo = struct('rMax', nan, 'rHatByParam', struct(), 'nThin', nThin, ...
  'keepChains', [], 'elapsedSec', nan, 'chains', [], ...
  'thetaRMax', nan, 'thetaRHatByParticipant', [], 'thetaKeepChains', [], ...
  'rhatParticipants', [], 'nCollected', 0, 'batch', 0);

engine = char(p.Results.engine);
monitorParams = cellstr(monitorParams(:));
modelsDir = fileparts(mfilename('fullpath'));
generalDir = fullfile(modelsDir, '..', 'general');
if strcmpi(engine, 'stan')
  figuresDir = fullfile(modelsDir, 'stan', 'figures');
  storageDir = fullfile(modelsDir, 'stan', 'storage');
  logDir = fullfile(modelsDir, 'stan', 'logs');
else
  figuresDir = fullfile(modelsDir, 'figures');
  storageDir = fullfile(modelsDir, 'storage');
  logDir = fullfile(modelsDir, 'logs');
end
if ~isfolder(logDir), mkdir(logDir); end
if ~isfolder(figuresDir), mkdir(figuresDir); end
checkpointFile = fullfile(logDir, sprintf('%s.checkpoint.mat', modelName));

addpath(generalDir);
cleanupObj = onCleanup(@() rmpath(generalDir)); %#ok<NASGU>

[data, d] = prepareIntertemporalChoiceData(p.Results.dataName, p.Results.dataDir);
if strcmpi(engine, 'stan')
  data.muPrecScale = p.Results.muPrecScale;
end
rhatParticipants = normalizeParticipantIndices( ...
  p.Results.rhatParticipants, double(d.nParticipants));
attemptInfo.rhatParticipants = rhatParticipants;
if isempty(rhatParticipants)
  fprintf('Convergence gate: all %d participants\n', double(d.nParticipants));
else
  fprintf('Convergence gate: participants [%s] (%d/%d); model still includes all participants\n', ...
    sprintf('%d ', rhatParticipants), numel(rhatParticipants), double(d.nParticipants));
end

if nargin < 3 || isempty(initGenerator)
  if any(strcmp(monitorParams, 'w'))
    initGenerator = @() struct('w', rand(d.nParticipants, 1) * 2 + 0.5);
  elseif any(strcmp(monitorParams, 'sigma'))
    initGenerator = @() struct('sigma', rand(d.nParticipants, 1) * 0.4 + 0.05);
  else
    initGenerator = @() struct('epsilon', rand(d.nParticipants, 1) * 0.4 + 0.05);
  end
end

fileName = sprintf('%s_%s_%s.mat', modelName, p.Results.dataName, engine);
storagePath = fullfile(storageDir, fileName);

if p.Results.figuresOnly
  if ~isfile(storagePath)
    error('runHierarchicalExecutionModel:missingStorage', ...
      'figuresOnly requested but no storage file: %s', storagePath);
  end
  S = load(storagePath, 'chains');
  showFinalHierarchicalFigures(S.chains, monitorParams, data, d, modelName, ...
    figuresDir, p.Results.saveFigures);
  converged = true;
  return;
end

if p.Results.preLoad && isfile(storagePath)
  fprintf('Checking stored chains: %s\n', storagePath);
  S = load(storagePath, 'chains', 'stats', 'diagnostics', 'info', 'rHatByParam');
  chains = S.chains;
  nStored = codaNSamples(chains);
  gateChains = tailCodaChains(chains, min(nStored, nSamplesMin));
  [rMax, keepChains, rHatByParam, rMaxAll] = hierarchicalSubsetRhat( ...
    gateChains, monitorParams, p.Results.keepChainsMin, p.Results.rhatCritical, ...
    rhatParticipants);
  attemptInfo.rMax = rMax;
  attemptInfo.rHatByParam = rHatByParam;
  attemptInfo.keepChains = keepChains;
  attemptInfo.nCollected = nStored;
  fprintf(['R-hat subset-max=%.3f (all-chain=%.3f) over {%s} (last %d); keep %d / %d; ' ...
    'stored n=%d (min %g)\n'], ...
    rMax, rMaxAll, strjoin(fieldnames(rHatByParam), ', '), codaNSamples(gateChains), ...
    numel(keepChains), p.Results.nChains, nStored, nSamplesMin);

  if rMax <= p.Results.rhatCritical && numel(keepChains) >= p.Results.keepChainsMin ...
      && nStored >= nSamplesMin
    fprintf('Stored chains meet R-hat <= %.3g and n >= %g — using saved fit\n', ...
      p.Results.rhatCritical, nSamplesMin);
    chains = subsetChainsHierarchical(chains, keepChains);
    showFinalHierarchicalFigures(chains, monitorParams, data, d, modelName, figuresDir, p.Results.saveFigures);
    converged = true;
    return;
  end
  fprintf('Stored chains do not meet gate (R-hat <= %.3g and n >= %g) — continuing\n', ...
    p.Results.rhatCritical, nSamplesMin);
end

if ~isfolder(storageDir)
  mkdir(storageDir);
end
if p.Results.resetThin || p.Results.coldStart
  clearMcmcState(storageDir, modelName);
  if isfile(checkpointFile)
    delete(checkpointFile);
    fprintf('Cleared sequential checkpoint %s\n', checkpointFile);
  end
end

[acc, lastInitCell, initFn, nBurnin, bestRhat, nCollected, batch] = ...
  resumeHierCheckpoint(checkpointFile, modelName, initGenerator, ...
    nBurninFirst, nBurninContinue);

% Seed from seedStoragePath or canonical storage (last checkpointKeep samples).
seedPath = char(p.Results.seedStoragePath);
if isempty(acc)
  seedFrom = '';
  if ~isempty(seedPath) && isfile(seedPath)
    seedFrom = seedPath;
  elseif ~p.Results.coldStart && isfile(storagePath)
    seedFrom = storagePath;
  end
  if ~isempty(seedFrom)
    fprintf('Seeding sequential state from: %s\n', seedFrom);
    Sseed = load(seedFrom, 'chains');
    acc = subsetChainsToN(Sseed.chains, p.Results.nChains);
    acc = trimCodaTail(acc, checkpointKeep);
    nCollected = codaNSamples(acc);
    lastInitCell = codaLastInits(acc, d.nParticipants, p.Results.nChains);
    initFn = makeSerialInitFn(lastInitCell);
    nBurnin = nBurninContinue;
    gateAcc = tailCodaChains(acc, min(nCollected, nSamplesMin));
    [bestRhat, ~, ~] = hierarchicalSubsetRhat( ...
      gateAcc, monitorParams, p.Results.keepChainsMin, p.Results.rhatCritical, ...
      rhatParticipants);
    batch = max(0, ceil(nCollected / nSamplesBatch));
    saveSeqCheckpoint(checkpointFile, modelName, acc, lastInitCell, ...
      bestRhat, nCollected, batch);
    fprintf('Seeded checkpoint: window n=%d batch=%d best trailing R-hat=%.3f\n', ...
      nCollected, batch, bestRhat);
  end
end

converged = false;
stats = struct();
diagnostics = struct();
info = struct();
rHatByParam = struct();
keepChains = [];
rMax = inf;
chains = [];

statusFile = fullfile(logDir, sprintf('%s.status', modelName));

while true
  batch = batch + 1;
  fprintf('\n--- %s | batch=%d | thin=%d | burnin=%g | collect=%g | stacked=%g ---\n', ...
    modelName, batch, nThin, nBurnin, nSamplesBatch, nCollected);
  fprintf('Wall time %s\n', datestr(now, 31));
  writeMcmcStatus(statusFile, sprintf( ...
    '%s | batch=%d | stacked=%g | collect +%g | best R-hat=%.3f | running', ...
    modelName, batch, nCollected, nSamplesBatch, bestRhat));

  [stats, chains, diagnostics, info] = fitHierarchicalExecutionOnce( ...
    engine, modelName, monitorParams, data, initFn, ...
    nThin, nBurnin, nSamplesBatch, p.Results);

  acc = appendCoda(acc, chains);
  acc = trimCodaTail(acc, checkpointKeep);
  if batch == 1 && nCollected == 0
    nCollected = codaNSamples(acc);
  else
    nCollected = nCollected + nSamplesBatch;
  end
  lastInitCell = codaLastInits(chains, d.nParticipants, p.Results.nChains);
  initFn = makeSerialInitFn(lastInitCell);
  nBurnin = nBurninContinue;

  gateAcc = tailCodaChains(acc, min(codaNSamples(acc), nSamplesMin));
  nGate = codaNSamples(gateAcc);
  [rMax, keepChains, rHatByParam, rMaxAll] = hierarchicalSubsetRhat( ...
    gateAcc, monitorParams, p.Results.keepChainsMin, p.Results.rhatCritical, ...
    rhatParticipants);
  fprintf(['stacked n=%d | trailing n=%d | subset R-hat max=%.3f (all-chain=%.3f) ' ...
    'over {%s}; keep %s (%d / %d)\n'], ...
    nCollected, nGate, rMax, rMaxAll, strjoin(fieldnames(rHatByParam), ', '), ...
    mat2str(keepChains), numel(keepChains), p.Results.nChains);

  attemptInfo.rMax = rMax;
  attemptInfo.rHatByParam = rHatByParam;
  attemptInfo.nThin = nThin;
  attemptInfo.keepChains = keepChains;
  attemptInfo.chains = acc;
  attemptInfo.nCollected = nCollected;
  attemptInfo.batch = batch;

  if p.Results.plotFailedRuns
    plotParticipantParameterCIs(gateAcc, monitorParams, ...
      'sgtitle', sprintf('%s | trailing=%d | max R-hat=%.3f', modelName, nGate, rMax), ...
      'savePath', hierarchicalFigurePath(figuresDir, modelName, ...
        sprintf('trail%d', nGate), p.Results.saveFigures));
  end

  try
    [thetaRMax, thetaRHatByP, thetaKeep] = assessHierarchicalThetaRhat( ...
      gateAcc, data, d, modelName, p.Results.keepChainsMin, ...
      p.Results.rhatCritical, rhatParticipants);
    attemptInfo.thetaRMax = thetaRMax;
    attemptInfo.thetaRHatByParticipant = thetaRHatByP;
    attemptInfo.thetaKeepChains = thetaKeep;
  catch ME
    fprintf('theta R-hat assessment skipped: %s\n', ME.message);
  end

  if rMax < bestRhat
    bestRhat = rMax;
    if p.Results.saveOnConverge
      saveHierBest(storagePath, gateAcc, keepChains, p.Results.keepChainsMin, ...
        stats, diagnostics, info, rHatByParam);
      fprintf('Updated best-so-far (trailing R-hat=%.3f, n=%d) → %s\n', ...
        bestRhat, nGate, storagePath);
    end
  end

  saveSeqCheckpoint(checkpointFile, modelName, acc, lastInitCell, ...
    bestRhat, nCollected, batch);

  converged = rMax <= p.Results.rhatCritical ...
    && numel(keepChains) >= p.Results.keepChainsMin ...
    && nCollected >= nSamplesMin;
  fprintf('Attempt done: converged=%d | trailing R-hat=%.3f | best=%.3f | stacked=%d (min %g)\n', ...
    converged, rMax, bestRhat, nCollected, nSamplesMin);
  writeMcmcStatus(statusFile, sprintf( ...
    '%s | batch=%d | stacked=%d | trail=%d | R-hat=%.3f | best=%.3f | keep=%d/%d | nMin=%g | converged=%d', ...
    modelName, batch, nCollected, nGate, rMax, bestRhat, numel(keepChains), ...
    p.Results.nChains, nSamplesMin, converged));

  if converged || p.Results.singleAttempt
    if ~converged
      fprintf('Not converged after batch %d (max R-hat=%.3f, n=%d)\n', ...
        batch, rMax, nCollected);
    end
    break;
  end
  fprintf('Appending another %g samples\n', nSamplesBatch);
end

if converged
  chains = subsetChainsHierarchical(gateAcc, keepChains);
  attemptInfo.chains = chains;
  if p.Results.saveOnConverge
    save(storagePath, 'chains', 'stats', 'diagnostics', 'info', 'rHatByParam', '-v7.3');
    fprintf('Saved %s (%d trailing samples)\n', storagePath, codaNSamples(chains));
    clearMcmcState(storageDir, modelName);
    if isfile(checkpointFile)
      delete(checkpointFile);
    end
    showFinalHierarchicalFigures(chains, monitorParams, data, d, modelName, ...
      figuresDir, p.Results.saveFigures);
  end
  grtable(chains, p.Results.rhatCritical);
  codatable(chains);
elseif p.Results.plotFailedRuns
  grtable(acc, p.Results.rhatCritical);
  codatable(acc);
end

end

function [stats, chains, diagnostics, info] = fitHierarchicalExecutionOnce( ...
    engine, modelName, monitorParams, data, initGenerator, ...
    nThin, nBurnin, nSamples, opts)

tic;
modelsDir = fileparts(mfilename('fullpath'));
if strcmpi(engine, 'stan')
  [modelPath, workDir] = stageStanEntropModel(modelsDir, modelName);
else
  modelPath = resolveJagsModelFile(modelsDir, ...
    sprintf('%s_%s.txt', modelName, engine));
  workDir = fullfile('tmp', modelName);
end
[stats, chains, diagnostics, info] = callbayes(engine, ...
  'model', modelPath, ...
  'data', data, ...
  'outputname', 'samples', ...
  'init', initGenerator, ...
  'datafilename', modelName, ...
  'initfilename', modelName, ...
  'scriptfilename', modelName, ...
  'logfilename', fullfile(workDir, 'log'), ...
  'nchains', opts.nChains, ...
  'nburnin', nBurnin, ...
  'nsamples', nSamples, ...
  'monitorparams', monitorParams, ...
  'thin', nThin, ...
  'workingdir', workDir, ...
  'verbosity', 0, ...
  'saveoutput', true, ...
  'allowunderscores', 1, ...
  'parallel', opts.doParallel);
fprintf('%s took %.1f s (thin=%d, burn-in=%g, collect=%g)\n', ...
  upper(engine), toc, nThin, nBurnin, nSamples);
end

function [acc, lastInitCell, initFn, nBurnin, bestRhat, nCollected, batch] = ...
    resumeHierCheckpoint(checkpointFile, modelName, meanInitFn, ...
    nBurninFirst, nBurninContinue)
acc = [];
lastInitCell = {};
initFn = meanInitFn;
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
  fprintf('Resuming from checkpoint: stacked=%d batch=%d best R-hat=%.3f\n', ...
    nCollected, batch, bestRhat);
end
end

function saveSeqCheckpoint(path, modelName, acc, lastInitCell, bestRhat, nCollected, batch)
save(path, 'modelName', 'acc', 'lastInitCell', 'bestRhat', 'nCollected', 'batch', '-v7.3');
end

function saveHierBest(storagePath, acc, keepChains, keepChainsMin, ...
    stats, diagnostics, info, rHatByParam)
if numel(keepChains) >= keepChainsMin
  chains = subsetChainsHierarchical(acc, keepChains);
else
  chains = acc;
end
save(storagePath, 'chains', 'stats', 'diagnostics', 'info', 'rHatByParam', '-v7.3');
end

function out = tailCodaChains(chains, nWin)
out = chains;
if isempty(chains) || nWin < 1
  return;
end
out = trimCodaTail(chains, nWin);
end

function out = trimCodaTail(chains, nKeep)
out = chains;
if isempty(chains) || nKeep < 1
  return;
end
fn = fieldnames(chains);
for k = 1:numel(fn)
  x = chains.(fn{k});
  if isnumeric(x) && size(x, 1) > nKeep
    out.(fn{k}) = x(end - nKeep + 1:end, :);
  end
end
end

function chainsOut = subsetChainsToN(chainsIn, nChains)
chainsOut = chainsIn;
if isempty(chainsIn) || nChains < 1
  return;
end
fn = fieldnames(chainsIn);
for k = 1:numel(fn)
  x = chainsIn.(fn{k});
  if isnumeric(x) && size(x, 2) > nChains
    chainsOut.(fn{k}) = x(:, 1:nChains);
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

function initCell = codaLastInits(chains, nP, nChains)
fn = fieldnames(chains);
initCell = cell(1, nChains);
scalarNames = {};
indexed = struct();
for k = 1:numel(fn)
  name = fn{k};
  if strcmp(name, 'deviance')
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

function indices = normalizeParticipantIndices(indices, nParticipants)
if isempty(indices)
  indices = [];
  return;
end
indices = unique(double(indices(:)'), 'stable');
if any(~isfinite(indices)) || any(indices ~= round(indices)) || ...
    any(indices < 1) || any(indices > nParticipants)
  error('runHierarchicalExecutionModel:badRhatParticipants', ...
    'rhatParticipants must contain integer indices from 1 to %d.', nParticipants);
end
end

function showFinalHierarchicalFigures(chains, monitorParams, data, d, modelName, figuresDir, saveFigures)
plotParticipantParameterCIs(chains, monitorParams, ...
  'sgtitle', sprintf('%s | final', modelName), ...
  'savePath', hierarchicalFigurePath(figuresDir, modelName, 'final', saveFigures));
plotHierarchicalPosteriorPredictiveFigure(chains, data, d, modelName, ...
  'sgtitle', sprintf('%s | posterior predictive | final', modelName), ...
  'savePath', hierarchicalFigurePath(figuresDir, modelName, 'postPredFinal', saveFigures));
printParameterMeans(chains, monitorParams);
end

function printParameterMeans(chains, monitorParams)
for k = 1:numel(monitorParams)
  pname = monitorParams{k};
  mu = get_matrix_from_coda(chains, pname, @mean);
  fprintf('%s mean (p1..p5): ', pname);
  fprintf('%.4g ', mu(1:min(5, numel(mu))));
  fprintf('\n');
end
end

function chainsOut = subsetChainsHierarchical(chainsIn, keepChains)
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

function [rMax, keepChains, rHatByParam, rMaxAll] = hierarchicalSubsetRhat( ...
  chains, monitorParams, keepChainsMin, rhatCritical, rhatParticipants)
% All-chain param breakdown for logs; gate statistic is subset max R-hat.
[rMaxAll, ~, rHatByParam] = maxRhatOverMonitoredParams( ...
  chains, monitorParams, keepChainsMin, rhatCritical, rhatParticipants);
[ok, keepChains, rMax] = findConvergedMonitoredChainSubset( ...
  chains, monitorParams, keepChainsMin, rhatCritical, rhatParticipants);
if ~ok && ~(isfinite(rMax))
  rMax = rMaxAll;
end
end
