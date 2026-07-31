%% runAlphaPriorReruns — re-fit SS/LL + latent-mixture entrop (faster schedule)
%
% Alpha is non-hierarchical (dbeta). Mixture hyperpriors use Trinity-safe
% camelCase names (muKappaEX, ...) so group means can be warm-started.
%
% Storage:
%   - Canonical mats archived once to *__alphaPriorOld.mat
%   - New fits write canonical paths used by drawFiguresEntrop
%   - storageTag = 'alphaPriorOld' loads archives in drawFiguresEntrop
%
% Mixture MCMC:
%   - Warm-start from hierarchical *ExecutionHierarchical_entrop + SS/LL
%   - thin 1..8; if still not converged, escalate burn-in then nSamples
%   - Best-so-far always written to canonical for drawFiguresEntrop
%
%   runAlphaPriorReruns

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
cd(modelsDir);
addpath(modelsDir);
addpath(fullfile(modelsDir, 'old'));

storageDir = fullfile(modelsDir, 'storage');
jagsDir = fullfile(modelsDir, 'jags');
jagsRobustDir = fullfile(modelsDir, 'jagsRobust');
if ~isfolder(storageDir), mkdir(storageDir); end
if ~isfolder(jagsRobustDir), mkdir(jagsRobustDir); end

dataName = 'intertemporalChoice';
engine = 'jags';
archiveTag = 'alphaPriorOld';
storageTagNew = '';
rhatCritical = 1.1;
keepChainsMin = 8;
nChains = 12;
nThinStart = 1;
maxThin = 8;
nBurninStart = 1e3;
maxBurnin = 8e3;
nSamplesStart = 5e3;
maxSamples = 2e4;
doSSLL = false;          % already re-fit after alpha prior change
doRefitUT = true;
doMixture = true;

baselineJags = fullfile(jagsDir, 'latentMixtureHierarchicalPrecision_entrop_jags.txt');
mixtureJobs = cell(3, 4);
mixtureJobs(1, :) = {'baseline', 1, ...
  'latentMixtureHierarchicalPrecision_entrop', baselineJags};
mixtureJobs(2, :) = {'half', 0.5, ...
  'latentMixtureHierarchicalMuPrecHalf_entrop', ...
  fullfile(jagsRobustDir, 'latentMixtureHierarchicalMuPrecHalf_entrop_jags.txt')};
mixtureJobs(3, :) = {'double', 2, ...
  'latentMixtureHierarchicalMuPrecDouble_entrop', ...
  fullfile(jagsRobustDir, 'latentMixtureHierarchicalMuPrecDouble_entrop_jags.txt')};

%% ---- archive current canonical results (once) ----
archiveModels = {'SS', 'LL', ...
  'unifiedTradeoffExecutionHierarchical_entrop', ...
  mixtureJobs{1, 3}, mixtureJobs{2, 3}, mixtureJobs{3, 3}};
for k = 1:numel(archiveModels)
  archiveStorageIfNeeded(storageDir, archiveModels{k}, dataName, engine, archiveTag);
end

%% ---- SS / LL ----
if doSSLL
  fprintf('\n=== Re-fitting SS / LL with dbeta(10,1) alpha ===\n');
  runContaminantModel('SS', dataName, engine, storageTagNew);
  runContaminantModel('LL', dataName, engine, storageTagNew);
end

%% ---- Re-fit hierarchical UT (standard priors, matching mixture) ----
if doRefitUT
  fprintf('\n=== Re-fitting unifiedTradeoffExecutionHierarchical_entrop (standard) ===\n');
  [~, dUt] = prepareIntertemporalChoiceData(dataName);
  initUt = @() struct('w', rand(dUt.nParticipants, 1) * 2 + 0.5);
  runHierarchicalExecutionModel( ...
    'unifiedTradeoffExecutionHierarchical_entrop', ...
    {'gamma', 'tau', 'kappa', 'vartheta', 'eta', 'w'}, ...
    initUt, ...
    'dataName', dataName, ...
    'preLoad', false, ...
    'rhatCritical', 1.4, ...
    'nBurnin', 2e3, ...
    'nSamples', 5e3, ...
    'nChains', 12, ...
    'saveFigures', true);
end

%% ---- latent mixture: baseline → half → double ----
if doMixture
  [~, d] = prepareIntertemporalChoiceData(dataName);
  zWarm = [];
  archiveMix = storageMatPath(storageDir, mixtureJobs{1, 3}, dataName, engine, archiveTag);
  if isfile(archiveMix)
    try
      S0 = load(archiveMix, 'chains');
      zWarm = round(get_matrix_from_coda(S0.chains, 'z', @mode));
      zWarm = max(1, min(11, zWarm(:)));
      fprintf('Warm-start z from archive %s\n', archiveMix);
    catch
    end
  end

  fprintf('\n=== Building mixture warm-start from hierarchical entrop fits ===\n');
  baseInitFn = makeLatentMixtureEntropInit(storageDir, dataName, engine, ...
    d.nParticipants, 'zMode', zWarm, 'jitterFrac', 0.05);

  for j = 1:size(mixtureJobs, 1)
    precLabel = mixtureJobs{j, 1};
    precScale = mixtureJobs{j, 2};
    modelName = mixtureJobs{j, 3};
    jagsPath = mixtureJobs{j, 4};

    fprintf('\n========== Mixture %s (%s) ==========\n', modelName, precLabel);

    if precScale ~= 1
      generateHierarchicalMuPrecJags(baselineJags, jagsPath, precScale);
    end

    canonicalPath = storageMatPath(storageDir, modelName, dataName, engine, storageTagNew);
    clearMcmcState(storageDir, modelName);

    nThin = nThinStart;
    nBurnin = nBurninStart;
    nSamples = nSamplesStart;
    bestRhat = inf;
    converged = false;
    initFn = baseInitFn;

    while true
      fprintf('\n--- %s | thin=%d | burnin=%g | samples=%g | full warm-start ---\n', ...
        modelName, nThin, nBurnin, nSamples);

      [converged, info] = runHierarchicalLatentMixtureModel(modelName, initFn, ...
        'dataName', dataName, ...
        'preLoad', false, ...
        'singleAttempt', true, ...
        'resetThin', true, ...
        'nThin', nThin, ...
        'nBurnin', nBurnin, ...
        'nSamples', nSamples, ...
        'nChains', nChains, ...
        'rhatCritical', rhatCritical, ...
        'keepChainsMin', keepChainsMin, ...
        'storageTag', storageTagNew, ...
        'saveOnConverge', true, ...
        'saveBestPath', canonicalPath, ...
        'previousBestRhat', bestRhat, ...
        'saveFigures', true);

      if isfield(info, 'bestRhat') && isfinite(info.bestRhat)
        bestRhat = info.bestRhat;
      end
      if isfield(info, 'chains') && ~isempty(info.chains)
        zWarm = round(get_matrix_from_coda(info.chains, 'z', @mode));
        zWarm = max(1, min(11, zWarm(:)));
        initFn = makeLatentMixtureEntropInit(storageDir, dataName, engine, ...
          d.nParticipants, 'zMode', zWarm, 'jitterFrac', 0.05);
      end

      fprintf('Attempt done: converged=%d | z R-hat=%.3f | best=%.3f | thin=%d burnin=%g samples=%g\n', ...
        converged, info.rMax, bestRhat, nThin, nBurnin, nSamples);

      if converged
        fprintf('=== %s reached z R-hat <= %.2f (thin=%d, burnin=%g, samples=%g) ===\n', ...
          modelName, rhatCritical, nThin, nBurnin, nSamples);
        break;
      end

      if nThin < maxThin
        nThin = nextThin(nThin);
        continue;
      end
      if nBurnin < maxBurnin
        nBurnin = min(maxBurnin, nBurnin * 2);
        fprintf('At maxThin=%d — escalating burn-in to %g\n', maxThin, nBurnin);
        continue;
      end
      if nSamples < maxSamples
        nSamples = min(maxSamples, nSamples * 2);
        fprintf('At max burn-in — escalating samples to %g\n', nSamples);
        continue;
      end

      fprintf('=== %s stopped (thin=%d, burnin=%g, samples=%g) without z R-hat <= %.2f (best=%.3f) ===\n', ...
        modelName, nThin, nBurnin, nSamples, rhatCritical, bestRhat);
      break;
    end
  end
end

fprintf('\n=== runAlphaPriorReruns finished ===\n');
fprintf('Canonical mats are current best for drawFiguresEntrop (storageTag = '''').\n');
fprintf('Previous results: storageTag = ''%s'' in drawFiguresEntrop.\n', archiveTag);

function runContaminantModel(modelName, dataName, engine, storageTag)
modelsDir = fileparts(mfilename('fullpath'));
storageDir = fullfile(modelsDir, 'storage');
generalDir = fullfile(modelsDir, '..', 'general');
addpath(generalDir);
cleanup = onCleanup(@() rmpath(generalDir));

[data, d] = prepareIntertemporalChoiceData(dataName);
data = struct('nParticipants', data.nParticipants, 'nTrials', data.nTrials, ...
  'decision', data.decision);

params = {'theta', 'alpha'};
nChains = 8;
nBurnin = 1e3;
nSamples = 1e3;
nThin = 1;
generator = @() struct('alpha', rand(d.nParticipants, 1) * 0.05 + 0.925);

storagePath = storageMatPath(storageDir, modelName, dataName, engine, storageTag);
modelPath = resolveJagsModelFile(modelsDir, sprintf('%s_%s.txt', modelName, engine));

fprintf('Fitting %s → %s\n', modelName, storagePath);
tic;
[stats, chains, diagnostics, info] = callbayes(engine, ...
  'model', modelPath, ...
  'data', data, ...
  'outputname', 'samples', ...
  'init', generator, ...
  'datafilename', modelName, ...
  'initfilename', modelName, ...
  'scriptfilename', modelName, ...
  'logfilename', fullfile('tmp', modelName), ...
  'nchains', nChains, ...
  'nburnin', nBurnin, ...
  'nsamples', nSamples, ...
  'monitorparams', params, ...
  'thin', nThin, ...
  'workingdir', fullfile('tmp', modelName), ...
  'verbosity', 0, ...
  'saveoutput', true, ...
  'allowunderscores', 1, ...
  'parallel', true);
fprintf('%s took %.1f s\n', upper(engine), toc);
save(storagePath, 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
fprintf('Saved %s\n', storagePath);
grtable(chains, 1.05);
codatable(chains);
end
