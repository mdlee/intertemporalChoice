%% runHierarchicalExecutionSequential — archive + sequential-batch refit
% Gate (best subset of at least 4 of 6 chains):
%   - thin = 1, batches of 1000
%   - max monitored R-hat < 1.05 on that subset
%   - at least 10,000 samples per chain
% Covers baseline + ITCH, then half/double prior-robustness variants.
% Hyperboloid (base + half/double) is scheduled last.
% Skips any mat that already meets the gate (preLoad); otherwise resumes
% checkpoint or seeds from existing storage.
%
% Live log: models/logs/runHierarchicalExecutionSequential.log
% Snapshot: models/logs/runHierarchicalExecutionSequential.status
% Stop file: models/logs/STOP_HIERARCHICAL_EXECUTION_SEQUENTIAL
%
%   runHierarchicalExecutionSequential

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
cd(modelsDir);
addpath(modelsDir);
addpath(fullfile(modelsDir, '..', 'general'));

logDir = fullfile(modelsDir, 'logs');
if ~isfolder(logDir), mkdir(logDir); end
diaryFile = fullfile(logDir, 'runHierarchicalExecutionSequential.log');
statusFile = fullfile(logDir, 'runHierarchicalExecutionSequential.status');
stopFile = fullfile(logDir, 'STOP_HIERARCHICAL_EXECUTION_SEQUENTIAL');
jagsRobustDir = fullfile(modelsDir, 'jagsRobust');
if ~isfolder(jagsRobustDir), mkdir(jagsRobustDir); end
if strcmpi(get(0, 'Diary'), 'on')
  diary off;
end
diary(diaryFile);
fprintf('Diary → %s\n', diaryFile);

storageDir = fullfile(modelsDir, 'storage');
dataName = 'intertemporalChoice';
engine = 'jags';
archiveTag = 'preSequentialBatch';
rhatCritical = 1.05;
nSamplesMin = 1e4;

% Hyperboloid last (base + half/double): thin=1 sequential is slowest / most brittle.
modelSpecs = {
  'exponentialExecutionHierarchical_entrop', ...
    {'kappa', 'w'};
  'hyperbolicExecutionHierarchical_entrop', ...
    {'kappa', 'w'};
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
  'hyperboloidExecutionHierarchical_entrop', ...
    {'kappa', 'tau', 'w'};
  };

precVariants = {
  '',      [];   % baseline
  'Half',  0.5;
  'Double', 2;
  };

% Build flat job list: baselines first, then half/double for each model.
jobs = {};
for m = 1:size(modelSpecs, 1)
  jobs(end+1, :) = {modelSpecs{m, 1}, modelSpecs{m, 2}, '', []}; %#ok<AGROW>
end
for m = 1:size(modelSpecs, 1)
  for pv = 2:size(precVariants, 1)
    jobs(end+1, :) = {modelSpecs{m, 1}, modelSpecs{m, 2}, ...
      precVariants{pv, 1}, precVariants{pv, 2}}; %#ok<AGROW>
  end
end
nJobs = size(jobs, 1);

writeHierStatus(statusFile, sprintf( ...
  ['starting | archiveTag=%s | %d jobs (8 baselines + half/double) | ' ...
   'R-hat<%.2f and n>=%g on >=4 of 6 chains'], ...
  archiveTag, nJobs, rhatCritical, nSamplesMin));

fprintf('=== Archive existing mats → storageTag=%s ===\n', archiveTag);
for j = 1:nJobs
  baseStem = jobs{j, 1};
  precLabel = jobs{j, 3};
  if isempty(precLabel)
    stem = baseStem;
  else
    stem = cognitiveMuPrecStem(baseStem, precLabel);
  end
  archiveStorageIfNeeded(storageDir, stem, dataName, engine, archiveTag);
end

[~, d] = prepareIntertemporalChoiceData(dataName);
initEntrop = @() struct('w', rand(d.nParticipants, 1) * 2 + 0.5);

for j = 1:nJobs
  if isfile(stopFile)
    msg = sprintf('stop file present — exiting before job %d/%d', j, nJobs);
    fprintf('%s\n', msg);
    writeHierStatus(statusFile, msg);
    diary off;
    return;
  end

  baseStem = jobs{j, 1};
  monitorParams = jobs{j, 2};
  precLabel = jobs{j, 3};
  precScale = jobs{j, 4};
  if isempty(precLabel)
    modelName = baseStem;
    jobLab = 'base';
  else
    modelName = cognitiveMuPrecStem(baseStem, precLabel);
    jobLab = sprintf('MuPrec%s', precLabel);
    srcJags = resolveJagsModelFile(modelsDir, sprintf('%s_jags.txt', baseStem));
    dstJags = fullfile(jagsRobustDir, sprintf('%s_jags.txt', modelName));
    generateHierarchicalMuPrecJags(srcJags, dstJags, precScale);
  end

  fprintf('\n========== [%d/%d] %s (%s) ==========\n', j, nJobs, modelName, jobLab);
  writeHierStatus(statusFile, sprintf( ...
    '[%d/%d] %s | %s | R-hat<%.2f and n>=%g on >=4/6 chains | fitting', ...
    j, nJobs, modelName, jobLab, rhatCritical, nSamplesMin));

  fitArgs = {
    'dataName', dataName, ...
    'preLoad', true, ...
    'resetThin', false, ...
    'keepChainsMin', 4, ...
    'rhatCritical', rhatCritical, ...
    'nSamplesMin', nSamplesMin, ...
    'nSamples', 1e3, ...
    'nChains', 6, ...
    'saveFigures', true};
  % Tradeoff already cold-started this overnight; resume checkpoint.
  % Gate is now best 4 of 6 chains (see findConvergedMonitoredChainSubset).
  % Hyperboloid: discard drifted sequential state; warm-start from thinned archive.
  if strcmp(modelName, 'hyperboloidExecutionHierarchical_entrop') && isempty(precLabel)
    seedPath = storageMatPath(storageDir, modelName, dataName, engine, archiveTag);
    fitArgs = [fitArgs, {
      'coldStart', true, ...
      'seedStoragePath', seedPath, ...
      'nBurnin', 2e3, ...
      'nBurninContinue', 5e2}];
    fprintf('Hyperboloid: cold restart, seed from %s\n', seedPath);
  end

  runHierarchicalExecutionModel( ...
    modelName, ...
    monitorParams, ...
    initEntrop, ...
    fitArgs{:});

  writeHierStatus(statusFile, sprintf('[%d/%d] %s | done', j, nJobs, modelName));
end

fprintf('\n=== runHierarchicalExecutionSequential finished %s ===\n', datestr(now, 31));
fprintf('Old chains: storageTag = ''%s''.\n', archiveTag);
fprintf('New chains: canonical (storageTag = '''').\n');
writeHierStatus(statusFile, sprintf('finished at %s', datestr(now, 31)));
diary off;

function writeHierStatus(statusFile, msg)
writeMcmcStatus(statusFile, msg);
end

function stem = cognitiveMuPrecStem(baseStem, precLabel)
stem = regexprep(baseStem, '_entrop$', ['MuPrec' precLabel '_entrop']);
if strcmp(stem, baseStem)
  error('cognitiveMuPrecStem:badStem', ...
    'Expected stem ending in _entrop, got: %s', baseStem);
end
end
