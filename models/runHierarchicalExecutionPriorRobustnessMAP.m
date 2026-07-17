%% runHierarchicalExecutionPriorRobustnessMAP
% Half/double mu_* prior fits for cognitive models that were MAP for at
% least one participant under the original mixture (Ex, UT, IT).
% All participants remain in every hierarchical fit, but convergence is
% assessed only for participants whose original-prior MAP is that model.
%
% Companion to runHierarchicalExecutionPriorRobustness.m for a second
% machine on the same Dropbox folder. Designed not to clash:
%   - Only MAP models (exponential already done → omitted here)
%   - Skips if storage/*.mat already exists
%   - Skips if *_mcmcState.mat exists (fit in progress elsewhere)
%   - Does not rewrite jagsRobust/*.txt if the file already exists
%
%   runHierarchicalExecutionPriorRobustnessMAP

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
jagsRobustDir = fullfile(modelsDir, 'jagsRobust');
storageDir = fullfile(modelsDir, 'storage');
if ~isfolder(jagsRobustDir), mkdir(jagsRobustDir); end
addpath(modelsDir);

dataName = 'intertemporalChoice';
engine = 'jags';
rhatCritical = 1.2;
keepChainsMin = 8;

% Original-mixture MAP cognitive models only (Ex already fitted elsewhere).
% Contaminants Gu/LL/SS have no mu_* hierarchical priors.
modelSpecs = {
  % 'exponentialExecutionHierarchical_entrop', ...
  %   {'kappa', 'w'}, 1;
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
initEntrop = @() struct('w', rand(d.nParticipants, 1) * 2 + 0.5);
mixturePath = fullfile(storageDir, ...
  sprintf('latentMixtureHierarchicalPrecision_entrop_%s_%s.mat', dataName, engine));
if ~isfile(mixturePath)
  error('Missing original-prior mixture fit: %s', mixturePath);
end
mixS = load(mixturePath, 'chains');
mapModel = originalMapModels(mixS.chains, double(d.nParticipants), 11);
clear mixS;

nRuns = size(modelSpecs, 1) * size(precVariants, 1);
fprintf('=== MAP-only cognitive prior robustness: %d candidate fits ===\n', nRuns);
fprintf('  models: %s\n', strjoin(modelSpecs(:, 1), ', '));
fprintf('  skip if storage or mcmcState already present (Dropbox-safe)\n');
for m = 1:size(modelSpecs, 1)
  pp = find(mapModel == modelSpecs{m, 3});
  fprintf('  %s convergence participants: [%s] (%d)\n', ...
    modelSpecs{m, 1}, sprintf('%d ', pp), numel(pp));
end

runIdx = 0;
nSkipped = 0;
nFitted = 0;
for m = 1:size(modelSpecs, 1)
  baseStem = modelSpecs{m, 1};
  monitorParams = modelSpecs{m, 2};
  rhatParticipants = find(mapModel == modelSpecs{m, 3});
  if isempty(rhatParticipants)
    fprintf('\n--- %s — skip (no original-prior MAP participants) ---\n', baseStem);
    nSkipped = nSkipped + size(precVariants, 1);
    runIdx = runIdx + size(precVariants, 1);
    continue;
  end
  srcJags = resolveJagsModelFile(modelsDir, sprintf('%s_jags.txt', baseStem));

  for pv = 1:size(precVariants, 1)
    precLabel = precVariants{pv, 1};
    precScale = precVariants{pv, 2};
    modelName = cognitiveMuPrecStem(baseStem, precLabel);
    runIdx = runIdx + 1;

    storagePath = fullfile(storageDir, ...
      sprintf('%s_%s_%s.mat', modelName, dataName, engine));
    mcmcStatePath = fullfile(storageDir, sprintf('%s_mcmcState.mat', modelName));

    if isfile(storagePath)
      fprintf('\n--- [%d/%d] %s — skip (storage exists) ---\n', ...
        runIdx, nRuns, modelName);
      nSkipped = nSkipped + 1;
      continue;
    end
    if isfile(mcmcStatePath)
      fprintf('\n--- [%d/%d] %s — skip (mcmcState present; in progress elsewhere?) ---\n', ...
        runIdx, nRuns, modelName);
      nSkipped = nSkipped + 1;
      continue;
    end

    dstJags = fullfile(jagsRobustDir, sprintf('%s_jags.txt', modelName));
    if ~isfile(dstJags)
      generateHierarchicalMuPrecJags(srcJags, dstJags, precScale);
    else
      fprintf('Using existing JAGS %s\n', dstJags);
    end

    fprintf('\n--- [%d/%d] %s | μ precision × %g ---\n', ...
      runIdx, nRuns, modelName, precScale);
    runHierarchicalExecutionModel(modelName, monitorParams, initEntrop, ...
      'dataName', dataName, ...
      'rhatCritical', rhatCritical, ...
      'keepChainsMin', keepChainsMin, ...
      'rhatParticipants', rhatParticipants, ...
      'preLoad', true);
    nFitted = nFitted + 1;
  end
end

fprintf('\n=== MAP-only prior-robustness pass complete ===\n');
fprintf('  fitted/attempted: %d | skipped: %d | storage: %s\n', ...
  nFitted, nSkipped, storageDir);

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
