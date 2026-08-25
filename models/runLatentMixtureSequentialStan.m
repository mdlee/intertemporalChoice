%% runLatentMixtureSequentialStan — 11-component mixture + μ-precision robustness (Stan)
% Jobs: baseline (muPrecScale=1), Half (0.5), Double (2). Same .stan file.
% Gate (best subset of at least 4 of 6 chains):
%   max participant z R-hat <= 1.10 on trailing 10,000 draws, n >= 10,000.
% Start after individual hierarchical Stan sequential finishes (see
% models/stan/logs/queueLatentMixtureAfterHierStan.sh).
%
% Live log: models/stan/logs/runLatentMixtureSequentialStan.log
% Snapshot: models/stan/logs/runLatentMixtureSequentialStan.status
% Stop file: models/stan/logs/STOP_LATENT_MIXTURE_SEQUENTIAL_STAN

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
cd(modelsDir);
addpath(modelsDir);
addpath(fullfile(modelsDir, '..', 'general'));

logDir = fullfile(modelsDir, 'stan', 'logs');
storageDir = fullfile(modelsDir, 'stan', 'storage');
if ~isfolder(logDir), mkdir(logDir); end
if ~isfolder(storageDir), mkdir(storageDir); end
diaryFile = fullfile(logDir, 'runLatentMixtureSequentialStan.log');
statusFile = fullfile(logDir, 'runLatentMixtureSequentialStan.status');
stopFile = fullfile(logDir, 'STOP_LATENT_MIXTURE_SEQUENTIAL_STAN');
if strcmpi(get(0, 'Diary'), 'on')
  diary off;
end
diary(diaryFile);
fprintf('Diary → %s\n', diaryFile);

dataName = 'intertemporalChoice';
archiveTag = 'preSequentialBatchStan';
rhatCritical = 1.10;
nSamplesMin = 1e4;

jobs = {
  'latentMixtureHierarchicalPrecision_entrop', 1, 'base';
  'latentMixtureHierarchicalPrecisionMuPrecHalf_entrop', 0.5, 'half';
  'latentMixtureHierarchicalPrecisionMuPrecDouble_entrop', 2, 'double';
  };
nJobs = size(jobs, 1);

baselineStems = {
  'exponentialExecutionHierarchical_entrop'
  'hyperbolicExecutionHierarchical_entrop'
  'hyperboloidExecutionHierarchical_entrop'
  'proportionalDifferencesExecutionHierarchical_entrop'
  'directDifferencesExecutionHierarchical_entrop'
  'tradeoffExecutionHierarchical_entrop'
  'unifiedTradeoffExecutionHierarchical_entrop'
  'itchExecutionHierarchical_entrop'
  };
for s = 1:numel(baselineStems)
  matPath = storageMatPath(storageDir, baselineStems{s}, dataName, 'stan', '');
  if ~isfile(matPath)
    error('runLatentMixtureSequentialStan:missingHier', ...
      ['Missing hierarchical Stan fit needed for mixture warm-start:\n  %s\n' ...
       'Wait until runHierarchicalExecutionSequentialStan finishes.'], matPath);
  end
end

writeMixStatus(statusFile, sprintf( ...
  ['starting | engine=stan | archiveTag=%s | %d jobs | ' ...
   'z R-hat<=%.2f and n>=%g on >=4 of 6 chains'], ...
  archiveTag, nJobs, rhatCritical, nSamplesMin));

fprintf('=== Archive existing Stan mixture mats → storageTag=%s ===\n', archiveTag);
for j = 1:nJobs
  archiveStorageIfNeeded(storageDir, jobs{j, 1}, dataName, 'stan', archiveTag);
end

for j = 1:nJobs
  if isfile(stopFile)
    msg = sprintf('stop file present — exiting before job %d/%d', j, nJobs);
    fprintf('%s\n', msg);
    writeMixStatus(statusFile, msg);
    diary off;
    return;
  end

  modelName = jobs{j, 1};
  precScale = jobs{j, 2};
  jobLab = jobs{j, 3};
  fprintf('\n========== [%d/%d] %s (%s) Stan mixture ==========\n', ...
    j, nJobs, modelName, jobLab);
  writeMixStatus(statusFile, sprintf( ...
    '[%d/%d] %s | %s | Stan mixture | z R-hat<=%.2f n>=%g on >=4/6 | fitting', ...
    j, nJobs, modelName, jobLab, rhatCritical, nSamplesMin));

  runHierarchicalLatentMixtureEntropStan(modelName, ...
    'dataName', dataName, ...
    'muPrecScale', precScale, ...
    'preLoad', true, ...
    'keepChainsMin', 4, ...
    'rhatCritical', rhatCritical, ...
    'nSamplesMin', nSamplesMin, ...
    'nSamples', 1e3, ...
    'nChains', 6, ...
    'saveFigures', true);

  writeMixStatus(statusFile, sprintf('[%d/%d] %s | Stan mixture | done', ...
    j, nJobs, modelName));
end

fprintf('\n=== runLatentMixtureSequentialStan finished %s ===\n', datestr(now, 31));
writeMixStatus(statusFile, sprintf('finished at %s', datestr(now, 31)));
diary off;

function writeMixStatus(statusFile, msg)
writeMcmcStatus(statusFile, msg);
end
