%% runLatentMixtureRobustnessSequential — half/double mixture, sequential batches
%
% Archives existing half/double mixture mats and JAGS, regenerates robustness
% JAGS from the current baseline (ITCH sign constraints; SS/LL dbeta(10,1);
% DD omega in (0,1)), then
% refits half then double the same way as runLatentMixtureSequential:
%   - thin = 1, 6 chains
%   - first batch: burn-in 1e3, collect 1e3
%   - later batches: last-state re-init, burn-in 200, collect 1e3, append z
%   - each batch redraws z from the regular mixture posterior (per person)
%   - R-hat gate 1.1: after each batch, if any subset of >=4 chains has
%     all participant z R-hats <= 1.1, save that subset and move on
%   - stop file checked at the start of each batch
%   - when both half and double finish, regenerate MS robustness figures
%
% Live log: models/logs/runLatentMixtureRobustnessSequential.log
% Snapshot: models/logs/runLatentMixtureRobustnessSequential.status
% Stop file: models/logs/STOP_LATENT_MIXTURE_ROBUSTNESS
% Checkpoint: models/logs/runLatentMixtureRobustnessSequential.checkpoint.mat

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
cd(modelsDir);
addpath(modelsDir);
addpath(fullfile(modelsDir, '..', 'general'));

logDir = fullfile(modelsDir, 'logs');
if ~isfolder(logDir), mkdir(logDir); end
diaryFile = fullfile(logDir, 'runLatentMixtureRobustnessSequential.log');
statusFile = fullfile(logDir, 'runLatentMixtureRobustnessSequential.status');
if strcmpi(get(0, 'Diary'), 'on')
  diary off;
end
diary(diaryFile);
fprintf('Diary → %s\n', diaryFile);

storageDir = fullfile(modelsDir, 'storage');
jagsDir = fullfile(modelsDir, 'jags');
jagsRobustDir = fullfile(modelsDir, 'jagsRobust');
jagsRobustOldDir = fullfile(jagsRobustDir, 'old');
figuresDir = fullfile(modelsDir, 'figures');
if ~isfolder(jagsRobustDir), mkdir(jagsRobustDir); end
if ~isfolder(jagsRobustOldDir), mkdir(jagsRobustOldDir); end

dataName = 'intertemporalChoice';
engine = 'jags';
archiveTag = 'preItchSign';
nModels = 11;

nChains = 6;
keepChainsMin = 4; % accept any converged subset of at least 4 chains
rhatCritical = 1.1;
nBurninFirst = 1e3;
nBurninContinue = 2e2;
nSamplesBatch = 1e3;
nThin = 1;
stopFile = fullfile(logDir, 'STOP_LATENT_MIXTURE_ROBUSTNESS');
checkpointFile = fullfile(logDir, 'runLatentMixtureRobustnessSequential.checkpoint.mat');

baselineName = 'latentMixtureHierarchicalPrecision_entrop';
baselineJags = fullfile(jagsDir, 'latentMixtureHierarchicalPrecision_entrop_jags.txt');
baselineMat = storageMatPath(storageDir, baselineName, dataName, engine, '');

mixtureJobs = {
  'half', 0.5, ...
    'latentMixtureHierarchicalMuPrecHalf_entrop', ...
    fullfile(jagsRobustDir, 'latentMixtureHierarchicalMuPrecHalf_entrop_jags.txt');
  'double', 2, ...
    'latentMixtureHierarchicalMuPrecDouble_entrop', ...
    fullfile(jagsRobustDir, 'latentMixtureHierarchicalMuPrecDouble_entrop_jags.txt')
  };

fprintf('=== Archive existing robustness mixture mats → storageTag=%s ===\n', archiveTag);
writeMixtureStatus(statusFile, sprintf('archiving robustness mats → %s', archiveTag));
for j = 1:size(mixtureJobs, 1)
  archiveStorageIfNeeded(storageDir, mixtureJobs{j, 3}, dataName, engine, archiveTag);
end

fprintf('=== Archive existing robustness mixture JAGS ===\n');
for j = 1:size(mixtureJobs, 1)
  srcJags = mixtureJobs{j, 4};
  [~, stem, ext] = fileparts(srcJags);
  dstJags = fullfile(jagsRobustOldDir, sprintf('%s_%s%s', stem, archiveTag, ext));
  if isfile(srcJags) && ~isfile(dstJags)
    copyfile(srcJags, dstJags);
    fprintf('Archived %s → %s\n', srcJags, dstJags);
  elseif isfile(dstJags)
    fprintf('JAGS archive already present: %s\n', dstJags);
  end
end

fprintf('=== Regenerate half/double JAGS from current baseline ===\n');
for j = 1:size(mixtureJobs, 1)
  generateHierarchicalMuPrecJags(baselineJags, mixtureJobs{j, 4}, mixtureJobs{j, 2});
end

[data, d] = prepareIntertemporalChoiceData(dataName);
data.nModels = nModels;

fprintf('\n=== Loading regular mixture z posterior for inits ===\n');
if ~isfile(baselineMat)
  error('runLatentMixtureRobustnessSequential:noBaseline', ...
    'Regular mixture mat missing: %s', baselineMat);
end
Zreg = loadMixtureZDraws(baselineMat, d.nParticipants);
fprintf('Regular z draws: %d samples × %d participants from %s\n', ...
  size(Zreg, 1), size(Zreg, 2), baselineMat);

fprintf('\n=== Building mean warm-start from hierarchical fits ===\n');
meanInitFn = makeLatentMixtureEntropInit(storageDir, dataName, engine, ...
  d.nParticipants, ...
  'zMode', [], ...
  'jitterFrac', 0.05, ...
  'summaryFn', @mean, ...
  'seed', 1);
meanInitFn = wrapInitWithSampledZ(meanInitFn, Zreg, nModels);

s0 = meanInitFn();
monitorParams = fieldnames(s0);
fprintf('Init fields (%d). z(1:5)=%s  omegaDD p1=%.4g  muOmegaDD=%.4g\n', ...
  numel(monitorParams), mat2str(s0.z(1:min(5,end))'), ...
  s0.omegaDD(1), s0.muOmegaDD);

for j = 1:size(mixtureJobs, 1)
  precLabel = mixtureJobs{j, 1};
  modelName = mixtureJobs{j, 3};

  fprintf('\n========== Mixture %s (%s) ==========\n', modelName, precLabel);
  writeMixtureStatus(statusFile, sprintf('starting %s (%s)', modelName, precLabel));

  canonicalPath = storageMatPath(storageDir, modelName, dataName, engine, '');
  clearMcmcState(storageDir, modelName);

  [zAcc, lastInitCell, initFn, nBurnin, bestRhat, nCollected, batch] = ...
    resumeRobustnessCheckpoint(checkpointFile, modelName, ...
      nChains, meanInitFn, nBurninFirst, nBurninContinue, rhatCritical);
  if ~isempty(lastInitCell)
    lastInitCell = overlaySampledZ(lastInitCell, Zreg, nModels);
    initFn = makeSerialInitFn(lastInitCell);
  end

  stopRequested = false;
  while true
    if isfile(stopFile)
      fprintf('\n=== Stop file present (%s) — ending after stacked=%d ===\n', ...
        stopFile, nCollected);
      stopRequested = true;
      break;
    end
    batch = batch + 1;
    fprintf('\n--- %s | batch=%d | thin=%d | burnin=%g | collect=%g | stacked=%g ---\n', ...
      modelName, batch, nThin, nBurnin, nSamplesBatch, nCollected);
    fprintf('Wall time %s\n', datestr(now, 31));
    writeMixtureStatus(statusFile, sprintf( ...
      '%s | batch=%d | stacked=%g | collect +%g | best z R-hat=%.3f | running', ...
      modelName, batch, nCollected, nSamplesBatch, bestRhat));

    tic;
    modelPath = resolveJagsModelFile(modelsDir, sprintf('%s_%s.txt', modelName, engine));
    [stats, chains, diagnostics, info] = callbayes(engine, ...
      'model', modelPath, ...
      'data', data, ...
      'outputname', 'samples', ...
      'init', initFn, ...
      'datafilename', modelName, ...
      'initfilename', modelName, ...
      'scriptfilename', modelName, ...
      'logfilename', fullfile('tmp', modelName), ...
      'nchains', nChains, ...
      'nburnin', nBurnin, ...
      'nsamples', nSamplesBatch, ...
      'monitorparams', monitorParams, ...
      'thin', nThin, ...
      'workingdir', fullfile('tmp', modelName), ...
      'verbosity', 0, ...
      'saveoutput', true, ...
      'allowunderscores', 1, ...
      'parallel', true);
    fprintf('%s took %.1f s (thin=%d, burn-in=%g, collect=%g)\n', ...
      upper(engine), toc, nThin, nBurnin, nSamplesBatch);

    zBatch = extractZCoda(chains);
    zAcc = appendCoda(zAcc, zBatch);
    nCollected = codaNSamples(zAcc);
    lastInitCell = overlaySampledZ( ...
      codaLastInits(chains, d.nParticipants, nChains, nModels), Zreg, nModels);
    initFn = makeSerialInitFn(lastInitCell);
    nBurnin = nBurninContinue;

    % All-chain diagnostics (for logging)
    [zRhatAll, ~, rHatEachAll] = maxRhatOverParticipants( ...
      zAcc, nChains, rhatCritical);
    fprintf('stacked n=%d | all-chain z R-hat: max=%.3f, min=%.3f, %d/%d at 1\n', ...
      nCollected, zRhatAll, min(rHatEachAll), sum(rHatEachAll <= 1 + 1e-9), numel(rHatEachAll));
    reportZRhatEach(rHatEachAll, rhatCritical);

    % Accept largest subset of size >= keepChainsMin with all R-hat <= gate
    [converged, keepChains, zRhat, rHatEach] = findConvergedZChainSubset( ...
      zAcc, keepChainsMin, rhatCritical);
    fprintf('subset gate (>=%d chains, R-hat<=%.2f): converged=%d | best max=%.3f | keep %s (%d/%d)\n', ...
      keepChainsMin, rhatCritical, converged, zRhat, mat2str(keepChains), ...
      numel(keepChains), nChains);
    worstLine = reportZRhatEach(rHatEach, rhatCritical);
    plotLatentMixtureZPosterior(zAcc, ...
      'nModels', nModels, ...
      'sgtitle', sprintf('%s | stacked=%d | max R-hat=%.3f', modelName, nCollected, zRhatAll), ...
      'savePath', hierarchicalFigurePath(figuresDir, modelName, ...
        sprintf('stacked%d', nCollected), true));

    if zRhat < bestRhat
      bestRhat = zRhat;
      saveMixtureBest(canonicalPath, zAcc, keepChains, keepChainsMin, ...
        stats, diagnostics, info);
      fprintf('Updated best-so-far (z R-hat=%.3f, n=%d, keep=%s) → %s\n', ...
        bestRhat, nCollected, mat2str(keepChains), canonicalPath);
    end

    fprintf('Attempt done: converged=%d | subset z R-hat=%.3f | best=%.3f | stacked=%d\n', ...
      converged, zRhat, bestRhat, nCollected);
    writeMixtureStatus(statusFile, sprintf( ...
      '%s | stacked=%d | subset z R-hat=%.3f | keep=%s | best=%.3f | converged=%d\n%s', ...
      modelName, nCollected, zRhat, mat2str(keepChains), bestRhat, converged, worstLine));
    saveSeqCheckpoint(checkpointFile, modelName, zAcc, lastInitCell, ...
      bestRhat, nCollected, batch);
    if converged
      saveMixtureBest(canonicalPath, zAcc, keepChains, keepChainsMin, ...
        stats, diagnostics, info);
      fprintf(['z subset of %d chains has R-hat <= %.2f (max=%.3f, keep=%s) — ' ...
        'saving and moving to the next robustness job\n'], ...
        numel(keepChains), rhatCritical, zRhat, mat2str(keepChains));
      break;
    end
    fprintf('Appending another %g samples\n', nSamplesBatch);
  end
  if stopRequested
    break;
  end
end

if stopRequested
  fprintf('\n=== runLatentMixtureRobustnessSequential PAUSED %s ===\n', datestr(now, 31));
  fprintf('Checkpoint left for resume. Remove STOP_LATENT_MIXTURE_ROBUSTNESS and re-run to continue.\n');
  writeMixtureStatus(statusFile, sprintf( ...
    'paused at %s | remove STOP_LATENT_MIXTURE_ROBUSTNESS and re-run to resume', ...
    datestr(now, 31)));
  diary off;
  return;
end

fprintf('\n=== runLatentMixtureRobustnessSequential finished %s ===\n', datestr(now, 31));
fprintf('Old chains: storageTag = ''%s''.\n', archiveTag);
fprintf('New chains: canonical (storageTag = '''').\n');
writeMixtureStatus(statusFile, sprintf('finished at %s', datestr(now, 31)));

% Regenerate MS robustness figures from the new half/double mats
fprintf('\n=== Regenerating robustness figures via drawFiguresEntrop ===\n');
writeMixtureStatus(statusFile, 'regenerating robustness figures');
resultsDir = fullfile(modelsDir, '..', 'results');
addpath(resultsDir);
storageTag = '';
analysisListOverride = {
  'latentMixturePosteriorsHalf';
  'latentMixturePosteriorsDouble';
  'latentMixturePriorRobustness';
  };
run(fullfile(resultsDir, 'drawFiguresEntrop.m'));
fprintf('=== Robustness figures done ===\n');
writeMixtureStatus(statusFile, sprintf('finished + figures at %s', datestr(now, 31)));
diary off;

function writeMixtureStatus(statusFile, msg)
writeMcmcStatus(statusFile, msg);
end

function fh = wrapInitWithSampledZ(baseFn, Z, nModels)
fh = @apply;
  function s = apply
    s = baseFn();
    s.z = sampleZFromMixturePosterior(Z, nModels);
  end
end

function initCell = overlaySampledZ(initCell, Z, nModels)
for c = 1:numel(initCell)
  initCell{c}.z = sampleZFromMixturePosterior(Z, nModels);
end
end

function zc = extractZCoda(chains)
fn = fieldnames(chains);
zc = struct();
for k = 1:numel(fn)
  if strncmp(fn{k}, 'z_', 2)
    zc.(fn{k}) = chains.(fn{k});
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
    if strcmp(base, 'z')
      v = max(1, min(nModels, round(v)));
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

function saveMixtureBest(canonicalPath, zAcc, keepChains, keepChainsMin, ...
    stats, diagnostics, info)
if numel(keepChains) >= keepChainsMin
  chains = subsetChainsByKeep(zAcc, keepChains);
else
  chains = zAcc;
end
save(canonicalPath, 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
end

function chainsOut = subsetChainsByKeep(chainsIn, keepChains)
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

function saveSeqCheckpoint(path, modelName, zAcc, lastInitCell, bestRhat, nCollected, batch)
save(path, 'modelName', 'zAcc', 'lastInitCell', 'bestRhat', 'nCollected', 'batch', '-v7.3');
end

function [zAcc, lastInitCell, initFn, nBurnin, bestRhat, nCollected, batch] = ...
    resumeRobustnessCheckpoint(checkpointFile, modelName, ...
      nChains, meanInitFn, nBurninFirst, nBurninContinue, rhatCritical)
zAcc = [];
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
  zAcc = ck.zAcc;
  lastInitCell = ck.lastInitCell;
  initFn = makeSerialInitFn(lastInitCell);
  nBurnin = nBurninContinue;
  bestRhat = ck.bestRhat;
  nCollected = ck.nCollected;
  batch = ck.batch;
  fprintf('Resuming from checkpoint: stacked=%d batch=%d best z R-hat=%.3f\n', ...
    nCollected, batch, bestRhat);
  [~, ~, rHatEach] = maxRhatOverParticipants(zAcc, nChains, rhatCritical);
  reportZRhatEach(rHatEach, rhatCritical);
end
end

function worstLine = reportZRhatEach(rHatEach, rhatCritical)
nP = numel(rHatEach);
[rSorted, ord] = sort(rHatEach(:), 'descend');
nOver = sum(rHatEach > rhatCritical);
fprintf('z R-hat by participant (worst first); %d/%d > %.2f:\n', ...
  nOver, nP, rhatCritical);
for k = 1:nP
  fprintf('  p%02d  %.3f\n', ord(k), rSorted(k));
end
nShow = min(5, nP);
parts = cell(1, nShow);
for k = 1:nShow
  parts{k} = sprintf('p%d=%.3f', ord(k), rSorted(k));
end
worstLine = sprintf('worst: %s | %d/%d > %.2f', ...
  strjoin(parts, ' '), nOver, nP, rhatCritical);
end
