%% runItchBetaSign — archive unconstrained ITCH MCMC and refit with sign constraints
% Sequential batches of 1000 until max R-hat <= 1.1 (runHierarchicalExecutionModel).
%
%   runItchBetaSign

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
cd(modelsDir);
addpath(modelsDir);
addpath(fullfile(modelsDir, '..', 'general'));

storageDir = fullfile(modelsDir, 'storage');
jagsDir = fullfile(modelsDir, 'jags');
jagsRobustDir = fullfile(modelsDir, 'jagsRobust');
dataName = 'intertemporalChoice';
engine = 'jags';
modelName = 'itchExecutionHierarchical_entrop';
archiveTag = 'preBetaSign';

archiveStorageIfNeeded(storageDir, modelName, dataName, engine, archiveTag);
archiveStorageIfNeeded(storageDir, 'itchExecutionHierarchicalMuPrecHalf_entrop', dataName, engine, archiveTag);
archiveStorageIfNeeded(storageDir, 'itchExecutionHierarchicalMuPrecDouble_entrop', dataName, engine, archiveTag);
clearMcmcState(storageDir, modelName);

fprintf('=== Regenerate ITCH half/double JAGS from constrained baseline ===\n');
srcJags = fullfile(jagsDir, 'itchExecutionHierarchical_entrop_jags.txt');
generateHierarchicalMuPrecJags(srcJags, ...
  fullfile(jagsRobustDir, 'itchExecutionHierarchicalMuPrecHalf_entrop_jags.txt'), 0.5);
generateHierarchicalMuPrecJags(srcJags, ...
  fullfile(jagsRobustDir, 'itchExecutionHierarchicalMuPrecDouble_entrop_jags.txt'), 2);

[~, d] = prepareIntertemporalChoiceData(dataName);
initFn = @() struct('w', rand(d.nParticipants, 1) * 2 + 0.5);

runHierarchicalExecutionModel( ...
  modelName, ...
  {'beta0', 'betaRA', 'betaRR', 'betaTA', 'betaTR', 'w'}, ...
  initFn, ...
  'dataName', dataName, ...
  'preLoad', false, ...
  'resetThin', true, ...
  'rhatCritical', 1.05, ...
  'nSamplesMin', 1e4, ...
  'saveFigures', true);

fprintf('\n=== runItchBetaSign finished %s ===\n', datestr(now, 31));
