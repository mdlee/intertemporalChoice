%% runHierarchicalExecutionPriorRobustness
% Fit each hierarchical cognitive model under half / double mu_* prior precision.
%
% Writes JAGS to models/jagsRobust/ and chains to models/storage/:
%   {base}MuPrecHalf_entrop_intertemporalChoice_jags.mat
%   {base}MuPrecDouble_entrop_intertemporalChoice_jags.mat
% where base is e.g. exponentialExecutionHierarchical (from ..._entrop).
%
% Used by results/drawFiguresEntrop → mapModelParameterRobustness.
%
%   runHierarchicalExecutionPriorRobustness

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
jagsDir = fullfile(modelsDir, 'jags');
jagsRobustDir = fullfile(modelsDir, 'jagsRobust');
storageDir = fullfile(modelsDir, 'storage');
if ~isfolder(jagsRobustDir), mkdir(jagsRobustDir); end
addpath(modelsDir);

dataName = 'intertemporalChoice';
engine = 'jags';
rhatCritical = 1.2;
keepChainsMin = 8;

% Same monitor lists as entrop hierarchical fits (epsilon → w).
modelSpecs = {
  'exponentialExecutionHierarchical_entrop', ...
    {'kappa', 'w'};
  'hyperbolicExecutionHierarchical_entrop', ...
    {'kappa', 'w'};
  'hyperboloidExecutionHierarchical_entrop', ...
    {'kappa', 'tau', 'w'};
  'proportionalDifferencesExecutionHierarchical_entrop', ...
    {'delta', 'w'};
  'directDifferencesExecutionHierarchical_entrop', ...
    {'delta', 'omega', 'w'};
  'tradeoffExecutionHierarchical_entrop', ...
    {'gamma', 'tau', 'kappa', 'vartheta', 'w'};
  'unifiedTradeoffExecutionHierarchical_entrop', ...
    {'gamma', 'tau', 'kappa', 'vartheta', 'eta', 'w'};
  'itchExecutionHierarchical_entrop', ...
    {'beta0', 'betaRA', 'betaRR', 'betaTA', 'betaTR', 'w'};
  };

precVariants = {
  'Half',  0.5;
  'Double', 2;
  };

[~, d] = prepareIntertemporalChoiceData(dataName);
initEntrop = @() struct('w', rand(d.nParticipants, 1) * 2 + 0.5);

nRuns = size(modelSpecs, 1) * size(precVariants, 1);
fprintf('=== Cognitive prior robustness: %d fits (%d models × %d precisions) ===\n', ...
  nRuns, size(modelSpecs, 1), size(precVariants, 1));

runIdx = 0;
for m = 1:size(modelSpecs, 1)
  baseStem = modelSpecs{m, 1};
  monitorParams = modelSpecs{m, 2};
  srcJags = resolveJagsModelFile(modelsDir, sprintf('%s_jags.txt', baseStem));

  for pv = 1:size(precVariants, 1)
    precLabel = precVariants{pv, 1};
    precScale = precVariants{pv, 2};
    modelName = cognitiveMuPrecStem(baseStem, precLabel);
    dstJags = fullfile(jagsRobustDir, sprintf('%s_jags.txt', modelName));

    generateHierarchicalMuPrecJags(srcJags, dstJags, precScale);

    runIdx = runIdx + 1;
    fprintf('\n--- [%d/%d] %s | μ precision × %g ---\n', ...
      runIdx, nRuns, modelName, precScale);
    runHierarchicalExecutionModel(modelName, monitorParams, initEntrop, ...
      'dataName', dataName, ...
      'rhatCritical', rhatCritical, ...
      'keepChainsMin', keepChainsMin, ...
      'preLoad', true);
  end
end

fprintf('\n=== Cognitive prior-robustness fits complete ===\n');
fprintf('Storage dir: %s\n', storageDir);

function stem = cognitiveMuPrecStem(baseStem, precLabel)
% cognitiveMuPrecStem  Insert MuPrecHalf/MuPrecDouble before _entrop.
stem = regexprep(baseStem, '_entrop$', ['MuPrec' precLabel '_entrop']);
if strcmp(stem, baseStem)
  error('cognitiveMuPrecStem:badStem', ...
    'Expected stem ending in _entrop, got: %s', baseStem);
end
end
