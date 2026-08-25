%% runHierarchicalExecutionSequentialStan — same sequential jobs, Trinity Stan
% Uses runHierarchicalExecutionModel(..., 'engine', 'stan').
% Gate (best subset of at least 4 of 6 chains):
%   - thin = 1, batches of 1000
%   - max monitored R-hat < 1.05 on that subset
%   - at least 10,000 samples per chain
% Half/double robustness: same .stan, data.muPrecScale = 0.5 or 2.
%
% Live log: models/stan/logs/runHierarchicalExecutionSequentialStan.log
% Snapshot: models/stan/logs/runHierarchicalExecutionSequentialStan.status
% Stop file: models/stan/logs/STOP_HIERARCHICAL_EXECUTION_SEQUENTIAL_STAN
%
%   runHierarchicalExecutionSequentialStan

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
cd(modelsDir);
addpath(modelsDir);
addpath(fullfile(modelsDir, '..', 'general'));

logDir = fullfile(modelsDir, 'stan', 'logs');
storageDir = fullfile(modelsDir, 'stan', 'storage');
if ~isfolder(logDir), mkdir(logDir); end
if ~isfolder(storageDir), mkdir(storageDir); end
diaryFile = fullfile(logDir, 'runHierarchicalExecutionSequentialStan.log');
statusFile = fullfile(logDir, 'runHierarchicalExecutionSequentialStan.status');
stopFile = fullfile(logDir, 'STOP_HIERARCHICAL_EXECUTION_SEQUENTIAL_STAN');
if strcmpi(get(0, 'Diary'), 'on')
  diary off;
end
diary(diaryFile);
fprintf('Diary → %s\n', diaryFile);

dataName = 'intertemporalChoice';
engine = 'stan';
archiveTag = 'preSequentialBatchStan';
rhatCritical = 1.05;
nSamplesMin = 1e4;

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
  '',      1;   % baseline
  'Half',  0.5;
  'Double', 2;
  };

jobs = {};
for m = 1:size(modelSpecs, 1)
  jobs(end+1, :) = {modelSpecs{m, 1}, modelSpecs{m, 2}, '', 1}; %#ok<AGROW>
end
for m = 1:size(modelSpecs, 1)
  for pv = 2:size(precVariants, 1)
    jobs(end+1, :) = {modelSpecs{m, 1}, modelSpecs{m, 2}, ...
      precVariants{pv, 1}, precVariants{pv, 2}}; %#ok<AGROW>
  end
end
nJobs = size(jobs, 1);

writeHierStatus(statusFile, sprintf( ...
  ['starting | engine=stan | archiveTag=%s | %d jobs | ' ...
   'R-hat<%.2f and n>=%g on >=4 of 6 chains'], ...
  archiveTag, nJobs, rhatCritical, nSamplesMin));

fprintf('=== Archive existing Stan mats → storageTag=%s ===\n', archiveTag);
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
  end

  fprintf('\n========== [%d/%d] %s (%s) Stan ==========\n', j, nJobs, modelName, jobLab);
  writeHierStatus(statusFile, sprintf( ...
    '[%d/%d] %s | %s | Stan | R-hat<%.2f and n>=%g on >=4/6 chains | fitting', ...
    j, nJobs, modelName, jobLab, rhatCritical, nSamplesMin));

  fitArgs = {
    'dataName', dataName, ...
    'engine', engine, ...
    'muPrecScale', precScale, ...
    'preLoad', true, ...
    'resetThin', false, ...
    'keepChainsMin', 4, ...
    'rhatCritical', rhatCritical, ...
    'nSamplesMin', nSamplesMin, ...
    'nSamples', 1e3, ...
    'nChains', 6, ...
    'saveFigures', true};

  runHierarchicalExecutionModel( ...
    modelName, ...
    monitorParams, ...
    initEntrop, ...
    fitArgs{:});

  writeHierStatus(statusFile, sprintf('[%d/%d] %s | Stan | done', j, nJobs, modelName));
end

fprintf('\n=== runHierarchicalExecutionSequentialStan finished %s ===\n', datestr(now, 31));
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
