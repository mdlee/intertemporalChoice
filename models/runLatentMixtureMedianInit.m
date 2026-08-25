%% runLatentMixtureMedianInit — baseline mixture at thin=1, median warm-start
%
% Initializes *all* mixture parameters from posterior medians of the
% hierarchical *ExecutionHierarchical_entrop (+ SS/LL) fits, then runs
% latentMixtureHierarchicalPrecision_entrop starting at thin=1.
%
% Usage (batch):
%   matlab -batch "run('runLatentMixtureMedianInit.m')"

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
cd(modelsDir);
addpath(modelsDir);
addpath(fullfile(modelsDir, '..', 'general'));

storageDir = fullfile(modelsDir, 'storage');
dataName = 'intertemporalChoice';
engine = 'jags';
modelName = 'latentMixtureHierarchicalPrecision_entrop';

rhatCritical = 1.1;
keepChainsMin = 8;
nChains = 12;
nThinStart = 1;
maxThin = 8;
nBurnin = 2e3;
nSamples = 5e3;

fprintf('=== %s | median warm-start | start thin=%d ===\n', modelName, nThinStart);
fprintf('Started %s\n', datestr(now, 31));

[~, d] = prepareIntertemporalChoiceData(dataName);

% Optional z warm-start from previous mixture (mode), else random in init.
zWarm = [];
prevMix = storageMatPath(storageDir, modelName, dataName, engine, '');
if isfile(prevMix)
  try
    S0 = load(prevMix, 'chains');
    zWarm = round(get_matrix_from_coda(S0.chains, 'z', @mode));
    zWarm = max(1, min(11, zWarm(:)));
    fprintf('Warm-start z from previous mixture %s\n', prevMix);
  catch ME
    fprintf('Could not load previous z (%s); using random z\n', ME.message);
  end
end

fprintf('\n=== Building init from individual-model posterior medians ===\n');
initFn = makeLatentMixtureEntropInit(storageDir, dataName, engine, ...
  d.nParticipants, ...
  'zMode', zWarm, ...
  'jitterFrac', 0.05, ...
  'summaryFn', @median, ...
  'seed', 1);

% Peek one draw so the log shows what we are initializing
s0 = initFn();
fn = fieldnames(s0);
fprintf('Init fields (%d): %s\n', numel(fn), strjoin(fn, ', '));
fprintf('  z(1:5)=%s ... z(end)=%d\n', mat2str(s0.z(1:min(5,end))'), s0.z(end));
fprintf('  w median-ish init p1..p5: ');
fprintf('%.3f ', s0.w(1:min(5,end)));
fprintf('\n');
if isfield(s0, 'muW')
  fprintf('  muW=%.3f  tauW=%.3f\n', s0.muW, s0.tauW);
end

canonicalPath = storageMatPath(storageDir, modelName, dataName, engine, '');
clearMcmcState(storageDir, modelName);

nThin = nThinStart;
bestRhat = inf;
converged = false;

while true
  fprintf('\n--- %s | thin=%d | burnin=%g | samples=%g | median warm-start ---\n', ...
    modelName, nThin, nBurnin, nSamples);
  fprintf('Wall time %s\n', datestr(now, 31));

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
    'storageTag', '', ...
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
      d.nParticipants, ...
      'zMode', zWarm, ...
      'jitterFrac', 0.05, ...
      'summaryFn', @median);
  end

  fprintf('Attempt done: converged=%d | z R-hat=%.3f | best=%.3f | thin=%d\n', ...
    converged, info.rMax, bestRhat, nThin);

  if converged
    fprintf('=== CONVERGED z R-hat <= %.2f (thin=%d) ===\n', rhatCritical, nThin);
    break;
  end
  if nThin >= maxThin
    fprintf('=== STOPPED at thin=%d without z R-hat <= %.2f (best=%.3f) ===\n', ...
      nThin, rhatCritical, bestRhat);
    break;
  end
  nThin = nextThin(nThin);
  fprintf('Not converged — next thin=%d\n', nThin);
end

fprintf('\n=== runLatentMixtureMedianInit finished %s ===\n', datestr(now, 31));
fprintf('Canonical path: %s\n', canonicalPath);
