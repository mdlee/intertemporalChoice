%% runHierarchicalExecutionPriorRobustnessMAP
% Half/double mu_* prior fits for cognitive models that were MAP for at
% least one participant under the original mixture (Ex, UT, IT).
% All participants remain in every hierarchical fit, but convergence is
% assessed only for participants whose original-prior MAP is that model.
%
% Draft-first schedule (for incremental mapModelParameterRobustness figures):
%   - Start at thin=1 with short MCMC (1e3 burn-in, 2e3 samples)
%   - Warm-start each chain near the converged *original-prior* hierarchical
%     posterior means for monitored participant parameters (5% jitter)
%   - After every attempt, save storage/*.mat if better (lower) R-hat than
%     any existing draft — even when not yet converged
%   - If a job does not converge, move on to the next incomplete job
%   - After one pass over incomplete jobs, double thin and repeat
%   - Stop when every job meets R-hat <= rhatCritical with enough chains
%
%   runHierarchicalExecutionPriorRobustnessMAP

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
resultsDir = fullfile(modelsDir, '..', 'results');
jagsRobustDir = fullfile(modelsDir, 'jagsRobust');
storageDir = fullfile(modelsDir, 'storage');
if ~isfolder(jagsRobustDir), mkdir(jagsRobustDir); end
addpath(modelsDir);
addpath(resultsDir);

dataName = 'intertemporalChoice';
engine = 'jags';
rhatCritical = 1.2;
keepChainsMin = 8;
nBurnin = 1e3;
nSamples = 2e3;
nThin = 1;

% Original-mixture MAP cognitive models (contaminants have no mu_* priors).
modelSpecs = {
  'exponentialExecutionHierarchical_entrop', ...
    {'kappa', 'w'}, 1;
  'unifiedTradeoffExecutionHierarchical_entrop', ...
    {'gamma', 'tau', 'kappa', 'vartheta', 'eta', 'w'}, 7;
  'itchExecutionHierarchical_entrop', ...
    {'beta0', 'betaRA', 'betaRR', 'betaTA', 'betaTR', 'w'}, 8;
  };

precVariants = {
  'Half',  0.5;
  'Double', 2;
  };

[~, d] = prepareIntertemporalChoiceData(dataName);
nParticipants = double(d.nParticipants);
mixturePath = fullfile(storageDir, ...
  sprintf('latentMixtureHierarchicalPrecision_entrop_%s_%s.mat', dataName, engine));
if ~isfile(mixturePath)
  error('Missing original-prior mixture fit: %s', mixturePath);
end
mixS = load(mixturePath, 'chains');
mapModel = originalMapModels(mixS.chains, nParticipants, 11);
clear mixS;

% Warm-start generators from converged original-prior hierarchical fits.
warmStartByBase = containers.Map('KeyType', 'char', 'ValueType', 'any');
for m = 1:size(modelSpecs, 1)
  baseStem = modelSpecs{m, 1};
  monitorParams = modelSpecs{m, 2};
  origPath = fullfile(storageDir, ...
    sprintf('%s_%s_%s.mat', baseStem, dataName, engine));
  if ~isfile(origPath)
    error('Missing original-prior hierarchical fit for warm start: %s', origPath);
  end
  fprintf('Warm-start source: %s\n', origPath);
  origS = load(origPath, 'chains');
  warmStartByBase(baseStem) = makeWarmStartInitFromChains( ...
    origS.chains, monitorParams, nParticipants, 'jitterFrac', 0.05);
  clear origS;
end

% Build job list: one entry per (MAP model × half/double) with ≥1 MAP participant.
jobs = buildMapJobs(modelSpecs, precVariants, mapModel, modelsDir, jagsRobustDir, ...
  storageDir, dataName, engine, warmStartByBase);
nJobs = numel(jobs);
if nJobs == 0
  error('No MAP cognitive jobs to run.');
end

fprintf('=== MAP prior-robustness draft loop: %d jobs ===\n', nJobs);
fprintf('  burn-in=%g | samples=%g | start thin=%d | R-hat<=%.3g | keepChains>=%d\n', ...
  nBurnin, nSamples, nThin, rhatCritical, keepChainsMin);
fprintf('  warm-start from original-prior hierarchical posterior means (5%% jitter)\n');
fprintf('  save every attempt if R-hat improves; advance thin after each full pass\n');
for j = 1:nJobs
  fprintf('  [%d] %s | MAP participants [%s] (%d)\n', ...
    j, jobs(j).modelName, sprintf('%d ', jobs(j).rhatParticipants), ...
    numel(jobs(j).rhatParticipants));
end

done = false(1, nJobs);
for j = 1:nJobs
  [rMax, keepChains] = storedFitRhat(jobs(j), keepChainsMin, rhatCritical);
  if isfinite(rMax) && rMax <= rhatCritical && numel(keepChains) >= keepChainsMin
    done(j) = true;
    fprintf('  already converged: %s (R-hat max=%.3f)\n', jobs(j).modelName, rMax);
  end
end

pass = 0;
while ~all(done)
  pass = pass + 1;
  incomplete = find(~done);
  fprintf('\n======== pass %d | thin=%d | incomplete %d/%d ========\n', ...
    pass, nThin, numel(incomplete), nJobs);

  improvedAny = false;
  for ii = 1:numel(incomplete)
    j = incomplete(ii);
    job = jobs(j);
    [oldRMax, ~] = storedFitRhat(job, keepChainsMin, rhatCritical);
    fprintf('\n--- [%d/%d incomplete] %s | thin=%d | stored R-hat=%s ---\n', ...
      ii, numel(incomplete), job.modelName, nThin, formatRhat(oldRMax));

    [converged, attemptInfo] = runHierarchicalExecutionModel( ...
      job.modelName, job.monitorParams, job.initGenerator, ...
      'dataName', dataName, ...
      'rhatCritical', rhatCritical, ...
      'keepChainsMin', keepChainsMin, ...
      'rhatParticipants', job.rhatParticipants, ...
      'nBurnin', nBurnin, ...
      'nSamples', nSamples, ...
      'nThin', nThin, ...
      'singleAttempt', true, ...
      'preLoad', false, ...
      'saveOnConverge', false, ...
      'plotFailedRuns', false, ...
      'saveFigures', false, ...
      'resetThin', true);

    newRMax = attemptInfo.rMax;
    keepChains = attemptInfo.keepChains;
    chains = attemptInfo.chains;
    rHatByParam = attemptInfo.rHatByParam;
    if isempty(chains)
      fprintf('No chains returned — leaving stored draft unchanged\n');
      continue;
    end

    isBetter = ~(isfinite(oldRMax)) || (isfinite(newRMax) && newRMax < oldRMax);
    meetGate = isfinite(newRMax) && newRMax <= rhatCritical && ...
      numel(keepChains) >= keepChainsMin;

    if isBetter || meetGate
      stats = struct();
      diagnostics = struct();
      info = struct('nThin', nThin, 'nBurnin', nBurnin, 'nSamples', nSamples, ...
        'rhatParticipants', job.rhatParticipants, 'draft', ~meetGate);
      if meetGate
        chains = subsetChainsHierarchicalLocal(chains, keepChains);
      end
      rMax = newRMax;
      save(job.storagePath, 'chains', 'stats', 'diagnostics', 'info', ...
        'rHatByParam', 'rMax', '-v7.3');
      improvedAny = true;
      if meetGate
        fprintf('Saved CONVERGED %s (R-hat max=%.3f, thin=%d)\n', ...
          job.storagePath, newRMax, nThin);
        clearMcmcState(storageDir, job.modelName);
        done(j) = true;
      else
        fprintf('Saved DRAFT %s (R-hat max=%.3f < stored %s; thin=%d)\n', ...
          job.storagePath, newRMax, formatRhat(oldRMax), nThin);
        saveUnsuccessfulMcmcState(storageDir, job.modelName, nThin);
      end
      refreshParameterRobustnessFigures(job.modelName, newRMax, nThin, meetGate);
    else
      fprintf('Kept stored draft (new R-hat=%.3f not better than %.3f)\n', ...
        newRMax, oldRMax);
      saveUnsuccessfulMcmcState(storageDir, job.modelName, nThin);
    end

    if converged
      done(j) = true;
    end
  end

  if all(done)
    break;
  end
  if ~improvedAny
    fprintf('No R-hat improvements this pass at thin=%d\n', nThin);
  end
  nThin = nextThin(nThin);
end

fprintf('\n=== MAP prior-robustness draft loop complete ===\n');
fprintf('  all %d jobs meet R-hat <= %.3g | storage: %s\n', nJobs, rhatCritical, storageDir);

function jobs = buildMapJobs(modelSpecs, precVariants, mapModel, modelsDir, ...
    jagsRobustDir, storageDir, dataName, engine, warmStartByBase)
jobs = struct('modelName', {}, 'monitorParams', {}, 'rhatParticipants', {}, ...
  'storagePath', {}, 'precScale', {}, 'initGenerator', {}, 'baseStem', {});
for m = 1:size(modelSpecs, 1)
  baseStem = modelSpecs{m, 1};
  monitorParams = modelSpecs{m, 2};
  rhatParticipants = find(mapModel == modelSpecs{m, 3});
  if isempty(rhatParticipants)
    fprintf('Skip %s — no original-prior MAP participants\n', baseStem);
    continue;
  end
  if ~isKey(warmStartByBase, baseStem)
    error('Missing warm-start generator for %s', baseStem);
  end
  srcJags = resolveJagsModelFile(modelsDir, sprintf('%s_jags.txt', baseStem));
  for pv = 1:size(precVariants, 1)
    precLabel = precVariants{pv, 1};
    precScale = precVariants{pv, 2};
    modelName = cognitiveMuPrecStem(baseStem, precLabel);
    dstJags = fullfile(jagsRobustDir, sprintf('%s_jags.txt', modelName));
    if ~isfile(dstJags)
      generateHierarchicalMuPrecJags(srcJags, dstJags, precScale);
    else
      fprintf('Using existing JAGS %s\n', dstJags);
    end
    job = struct();
    job.modelName = modelName;
    job.monitorParams = monitorParams;
    job.rhatParticipants = rhatParticipants;
    job.storagePath = fullfile(storageDir, ...
      sprintf('%s_%s_%s.mat', modelName, dataName, engine));
    job.precScale = precScale;
    job.initGenerator = warmStartByBase(baseStem);
    job.baseStem = baseStem;
    jobs(end + 1) = job; %#ok<AGROW>
  end
end
end

function [rMax, keepChains, rHatByParam] = storedFitRhat(job, keepChainsMin, rhatCritical)
rMax = inf;
keepChains = [];
rHatByParam = struct();
if ~isfile(job.storagePath)
  return;
end
S = load(job.storagePath, 'chains');
if ~isfield(S, 'chains') || isempty(S.chains)
  return;
end
[rMax, keepChains, rHatByParam] = maxRhatOverMonitoredParams( ...
  S.chains, job.monitorParams, keepChainsMin, rhatCritical, job.rhatParticipants);
end

function s = formatRhat(rMax)
if ~(isfinite(rMax))
  s = 'none';
else
  s = sprintf('%.3f', rMax);
end
end

function mapModel = originalMapModels(chains, nParticipants, nModels)
mapModel = nan(1, nParticipants);
zMatrices = codaIndexedMatrices(chains, 'z');
if numel(zMatrices) < nParticipants
  error('Expected z draws for %d participants, found %d.', ...
    nParticipants, numel(zMatrices));
end
for pp = 1:nParticipants
  z = zMatrices{pp}(:);
  probabilities = zeros(nModels, 1);
  for mi = 1:nModels
    probabilities(mi) = mean(z == mi);
  end
  [~, mapModel(pp)] = max(probabilities);
end
end

function stem = cognitiveMuPrecStem(baseStem, precLabel)
stem = regexprep(baseStem, '_entrop$', ['MuPrec' precLabel '_entrop']);
if strcmp(stem, baseStem)
  error('cognitiveMuPrecStem:badStem', ...
    'Expected stem ending in _entrop, got: %s', baseStem);
end
end

function chainsOut = subsetChainsHierarchicalLocal(chainsIn, keepChains)
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

function refreshParameterRobustnessFigures(modelName, rMax, nThin, converged)
% Rebuild git-linked repo eps/png and print a loud console note.
try
  fprintf('\n');
  fprintf('**************************************************************\n');
  fprintf('* PRIOR-ROBUSTNESS FIGURE REFRESH\n');
  fprintf('* trigger: improved fit for %s\n', modelName);
  fprintf('* R-hat max=%.3f | thin=%d | %s\n', ...
    rMax, nThin, ternaryLabel(converged, 'CONVERGED', 'DRAFT'));
  fprintf('**************************************************************\n');
  pathsParam = refreshMapModelParameterRobustnessFigures();
  pathsMatch = refreshMapModelForcedChoiceMatchFigures();
  fprintf('* Updated repo figures (git-linked from manuscript):\n');
  fprintf('*   %s\n', pathsParam.eps);
  fprintf('*   %s\n', pathsParam.png);
  fprintf('*   %s\n', pathsMatch.eps);
  fprintf('*   %s\n', pathsMatch.png);
  fprintf('**************************************************************\n\n');
catch ME
  fprintf(2, 'Prior-robustness figure refresh failed: %s\n', ME.message);
  if ~isempty(ME.stack)
    fprintf(2, '  at %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
  end
end
end

function s = ternaryLabel(tf, a, b)
if tf
  s = a;
else
  s = b;
end
end
