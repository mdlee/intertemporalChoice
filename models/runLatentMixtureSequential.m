%% runLatentMixtureSequential — mixture MCMC by appending 1000-draw batches
%
% Archives current canonical mixture mats to *__preItchSign.mat, then refits
% the baseline mixture with:
%   - thin = 1 always, 6 chains
%   - first batch: burn-in 1e3, collect 1e3
%   - later batches: re-init from last draw of each chain, short burn-in
%     (JAGS recompile/adapt), collect 1e3 more, APPEND to accumulated z
%   - R-hat is computed on the stacked samples (logged; does not stop)
%   - no sample cap; keep appending until models/logs/STOP_LATENT_MIXTURE exists
%     (checked at the start of each batch)
%   - warm-start first batch from posterior means of hierarchical fits
%   - each batch writes a checkpoint so a restart can resume the stack
%
% Half/double robustness mixture: runLatentMixtureRobustnessSequential.m
%
% Trinity starts a new JAGS process each batch, so this is continuation by
% last-state init + append, not a single JAGS update() stream.
%
% Live log: models/logs/runLatentMixtureSequential.log
% Snapshot: models/logs/runLatentMixtureSequential.status
% Stop file: models/logs/STOP_LATENT_MIXTURE
% Checkpoint: models/logs/runLatentMixtureSequential.checkpoint.mat

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
cd(modelsDir);
addpath(modelsDir);
addpath(fullfile(modelsDir, '..', 'general'));

logDir = fullfile(modelsDir, 'logs');
if ~isfolder(logDir), mkdir(logDir); end
diaryFile = fullfile(logDir, 'runLatentMixtureSequential.log');
statusFile = fullfile(logDir, 'runLatentMixtureSequential.status');
if strcmpi(get(0, 'Diary'), 'on')
  diary off;
end
diary(diaryFile);
fprintf('Diary → %s\n', diaryFile);

storageDir = fullfile(modelsDir, 'storage');
jagsDir = fullfile(modelsDir, 'jags');
jagsRobustDir = fullfile(modelsDir, 'jagsRobust');
figuresDir = fullfile(modelsDir, 'figures');
if ~isfolder(jagsRobustDir), mkdir(jagsRobustDir); end

dataName = 'intertemporalChoice';
engine = 'jags';
archiveTag = 'preItchSign';
nModels = 11;

nChains = 6;
keepChainsMin = 6;
rhatCritical = 1.1;
nBurninFirst = 1e3;
nBurninContinue = 2e2;
nSamplesBatch = 1e3;
nThin = 1;
stopFile = fullfile(logDir, 'STOP_LATENT_MIXTURE');
checkpointFile = fullfile(logDir, 'runLatentMixtureSequential.checkpoint.mat');

baselineJags = fullfile(jagsDir, 'latentMixtureHierarchicalPrecision_entrop_jags.txt');
% Baseline only. Half/double robustness: runLatentMixtureRobustnessSequential.m
mixtureJobs = {
  'baseline', 1, ...
    'latentMixtureHierarchicalPrecision_entrop', baselineJags;
  };

fprintf('=== Archive current mixture chains → storageTag=%s ===\n', archiveTag);
writeMixtureStatus(statusFile, sprintf('archiving canonical mixture mats → %s', archiveTag));
for j = 1:size(mixtureJobs, 1)
  archiveStorageIfNeeded(storageDir, mixtureJobs{j, 3}, dataName, engine, archiveTag);
end

fprintf('=== Regenerate half/double JAGS from current baseline ===\n');
generateHierarchicalMuPrecJags(baselineJags, ...
  fullfile(jagsRobustDir, 'latentMixtureHierarchicalMuPrecHalf_entrop_jags.txt'), 0.5);
generateHierarchicalMuPrecJags(baselineJags, ...
  fullfile(jagsRobustDir, 'latentMixtureHierarchicalMuPrecDouble_entrop_jags.txt'), 2);

[data, d] = prepareIntertemporalChoiceData(dataName);
data.nModels = nModels;

fprintf('\n=== Building mean warm-start from hierarchical fits ===\n');
meanInitFn = makeLatentMixtureEntropInit(storageDir, dataName, engine, ...
  d.nParticipants, ...
  'zMode', [], ...
  'jitterFrac', 0.05, ...
  'summaryFn', @mean, ...
  'seed', 1);

s0 = meanInitFn();
monitorParams = fieldnames(s0);
fprintf('Init fields (%d). z(1:5)=%s  omegaDD p1=%.4g  muOmegaDD=%.4g\n', ...
  numel(monitorParams), mat2str(s0.z(1:min(5,end))'), ...
  s0.omegaDD(1), s0.muOmegaDD);

for j = 1:size(mixtureJobs, 1)
  precLabel = mixtureJobs{j, 1};
  precScale = mixtureJobs{j, 2};
  modelName = mixtureJobs{j, 3};
  jagsPath = mixtureJobs{j, 4};

  fprintf('\n========== Mixture %s (%s) ==========\n', modelName, precLabel);
  writeMixtureStatus(statusFile, sprintf('starting %s (%s)', modelName, precLabel));
  if precScale ~= 1
    generateHierarchicalMuPrecJags(baselineJags, jagsPath, precScale);
  end

  canonicalPath = storageMatPath(storageDir, modelName, dataName, engine, '');
  clearMcmcState(storageDir, modelName);

  [zAcc, lastInitCell, initFn, nBurnin, bestRhat, nCollected, batch] = ...
    resumeMixtureCheckpoint(checkpointFile, modelName, canonicalPath, ...
      fullfile(modelsDir, 'tmp', modelName), nChains, meanInitFn, ...
      nBurninFirst, nBurninContinue, rhatCritical);

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
    lastInitCell = codaLastInits(chains, d.nParticipants, nChains, nModels);
    initFn = makeSerialInitFn(lastInitCell);
    nBurnin = nBurninContinue;

    [zRhat, keepChains, rHatEach] = maxRhatOverParticipants( ...
      zAcc, keepChainsMin, rhatCritical);
    fprintf('stacked n=%d | z R-hat: max=%.3f, min=%.3f, %d/%d at 1; keep %d / %d chains\n', ...
      nCollected, zRhat, min(rHatEach), sum(rHatEach <= 1 + 1e-9), numel(rHatEach), ...
      numel(keepChains), nChains);
    worstLine = reportZRhatEach(rHatEach, rhatCritical);
    try
      plotLatentMixtureZPosterior(zAcc, ...
        'nModels', nModels, ...
        'sgtitle', sprintf('%s | stacked=%d | max R-hat=%.3f', modelName, nCollected, zRhat), ...
        'savePath', hierarchicalFigurePath(figuresDir, modelName, ...
          sprintf('stacked%d', nCollected), true));
    catch ME
      fprintf('Figure save skipped: %s\n', ME.message);
    end

    if zRhat < bestRhat
      bestRhat = zRhat;
      saveMixtureBest(canonicalPath, zAcc, keepChains, keepChainsMin, ...
        stats, diagnostics, info);
      fprintf('Updated best-so-far (z R-hat=%.3f, n=%d) → %s\n', ...
        bestRhat, nCollected, canonicalPath);
    end

    converged = zRhat <= rhatCritical && numel(keepChains) >= keepChainsMin;
    fprintf('Attempt done: converged=%d | z R-hat=%.3f | best=%.3f | stacked=%d\n', ...
      converged, zRhat, bestRhat, nCollected);
    writeMixtureStatus(statusFile, sprintf( ...
      '%s | stacked=%d | z R-hat=%.3f | best=%.3f | converged=%d\n%s', ...
      modelName, nCollected, zRhat, bestRhat, converged, worstLine));
    saveSeqCheckpoint(checkpointFile, modelName, zAcc, lastInitCell, ...
      bestRhat, nCollected, batch);
    if converged
      fprintf('z R-hat <= %.2f — keeping samples and continuing until stop file\n', ...
        rhatCritical);
    end
    fprintf('Appending another %g samples\n', nSamplesBatch);
  end
  if stopRequested
    break;
  end
end

fprintf('\n=== runLatentMixtureSequential finished %s ===\n', datestr(now, 31));
fprintf('Old chains: storageTag = ''%s''.\n', archiveTag);
fprintf('New chains: canonical (storageTag = '''').\n');
writeMixtureStatus(statusFile, sprintf('finished at %s', datestr(now, 31)));
diary off;

function writeMixtureStatus(statusFile, msg)
writeMcmcStatus(statusFile, msg);
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
    resumeMixtureCheckpoint(checkpointFile, modelName, canonicalPath, tmpDir, ...
      nChains, meanInitFn, nBurninFirst, nBurninContinue, rhatCritical)
zAcc = [];
lastInitCell = {};
initFn = meanInitFn;
nBurnin = nBurninFirst;
bestRhat = inf;
nCollected = 0;
batch = 0;

if isfile(checkpointFile)
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
    [~, ~, rHatEach] = maxRhatOverParticipants(zAcc, 6, rhatCritical);
    reportZRhatEach(rHatEach, rhatCritical);
    return;
  end
end

initCell = cell(1, nChains);
for c = 1:nChains
  p = fullfile(tmpDir, sprintf('%s_%d.init', modelName, c));
  if ~isfile(p)
    return;
  end
  try
    initCell{c} = parseJagsInitFile(p);
  catch
    return;
  end
end
if ~isfile(canonicalPath)
  return;
end
S = load(canonicalPath, 'chains');
zTry = extractZCoda(S.chains);
nZ = codaNSamples(zTry);
if nZ <= 0
  return;
end
zAcc = zTry;
lastInitCell = initCell;
initFn = makeSerialInitFn(lastInitCell);
nBurnin = nBurninContinue;
nCollected = nZ;
batch = round(nZ / 1e3);
fprintf('Resuming from canonical z (n=%d) + tmp last-state inits\n', nCollected);
[~, ~, rHatEach] = maxRhatOverParticipants(zAcc, 6, rhatCritical);
reportZRhatEach(rHatEach, rhatCritical);
end

function s = parseJagsInitFile(path)
txt = fileread(path);
txt = regexprep(txt, '#[^\n]*', '');
s = struct();
[starts, ends, tokens] = regexp(txt, '"([A-Za-z]\w*)"\s*<-\s*', ...
  'start', 'end', 'tokens');
for i = 1:numel(starts)
  name = tokens{i}{1};
  if i < numel(starts)
    chunk = strtrim(txt(ends(i)+1:starts(i+1)-1));
  else
    chunk = strtrim(txt(ends(i)+1:end));
  end
  s.(name) = parseJagsValue(chunk);
end
if isfield(s, 'z')
  s.z = max(1, round(s.z(:)));
end
end

function v = parseJagsValue(chunk)
chunk = strtrim(chunk);
if strncmp(chunk, 'c(', 2)
  inner = chunk(3:end);
  inner = regexprep(inner, '\)\s*$', '');
  v = sscanf(strrep(inner, ',', ' '), '%f');
  v = v(:);
else
  v = sscanf(chunk, '%f', 1);
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
