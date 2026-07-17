function [converged, attemptInfo] = runHierarchicalExecutionModel(modelName, monitorParams, initGenerator, varargin)
%RUNHIERARCHICALEXECUTIONMODEL  Fit one hierarchical execution model (all participants).
%
%   MCMC defaults: 12 chains (≥8 required for convergence),
%   burn-in 2e3 (fixed), 5e3 samples (fixed), thin 1 (×2 per failure).
%   Unless singleAttempt is true, retries with doubled thinning until converged.
%   preLoad (default true): use storage only if max R-hat <= rhatCritical.
%   Failed runs save storage/{modelName}_mcmcState.mat.
%
%   [CONVERGED, ATTEMPTINFO] = ... also returns rMax, rHatByParam, nThin, keepChains.

p = inputParser;
addParameter(p, 'preLoad', true, @islogical);
addParameter(p, 'singleAttempt', false, @islogical);
addParameter(p, 'figuresOnly', false, @islogical);
addParameter(p, 'plotFailedRuns', true, @islogical);
addParameter(p, 'rhatCritical', 1.2, @isnumeric);
addParameter(p, 'keepChainsMin', 8, @isnumeric);
addParameter(p, 'nChains', 12, @isnumeric);
addParameter(p, 'nBurnin', 2e3, @isnumeric);
addParameter(p, 'nSamples', 5e3, @isnumeric);
addParameter(p, 'nThin', 1, @isnumeric);
addParameter(p, 'doParallel', true, @islogical);
addParameter(p, 'dataName', 'intertemporalChoice', @ischar);
addParameter(p, 'dataDir', '', @ischar);
addParameter(p, 'saveFigures', true, @islogical);
addParameter(p, 'resetThin', false, @islogical);
addParameter(p, 'saveOnConverge', true, @islogical);
parse(p, varargin{:});

attemptInfo = struct('rMax', nan, 'rHatByParam', struct(), 'nThin', p.Results.nThin, ...
  'keepChains', [], 'elapsedSec', nan, 'chains', [], ...
  'thetaRMax', nan, 'thetaRHatByParticipant', [], 'thetaKeepChains', []);

engine = 'jags';
monitorParams = cellstr(monitorParams(:));
modelsDir = fileparts(mfilename('fullpath'));
generalDir = fullfile(modelsDir, '..', 'general');
figuresDir = fullfile(modelsDir, 'figures');
storageDir = fullfile(modelsDir, 'storage');

addpath(generalDir);
cleanupObj = onCleanup(@() rmpath(generalDir));

[data, d] = prepareIntertemporalChoiceData(p.Results.dataName, p.Results.dataDir);

if nargin < 3 || isempty(initGenerator)
  if any(strcmp(monitorParams, 'w'))
    initGenerator = @() struct('w', rand(d.nParticipants, 1) * 2 + 0.5);
  elseif any(strcmp(monitorParams, 'sigma'))
    initGenerator = @() struct('sigma', rand(d.nParticipants, 1) * 0.4 + 0.05);
  else
    initGenerator = @() struct('epsilon', rand(d.nParticipants, 1) * 0.4 + 0.05);
  end
end

fileName = sprintf('%s_%s_%s.mat', modelName, p.Results.dataName, engine);
storagePath = fullfile(storageDir, fileName);

if p.Results.figuresOnly
  if ~isfile(storagePath)
    error('runHierarchicalExecutionModel:missingStorage', ...
      'figuresOnly requested but no storage file: %s', storagePath);
  end
  S = load(storagePath, 'chains');
  showFinalHierarchicalFigures(S.chains, monitorParams, data, d, modelName, ...
    figuresDir, p.Results.saveFigures);
  converged = true;
  return;
end

if p.Results.preLoad && isfile(storagePath)
  fprintf('Checking stored chains: %s\n', storagePath);
  S = load(storagePath, 'chains', 'stats', 'diagnostics', 'info', 'rHatByParam');
  chains = S.chains;
  if isfield(S, 'stats'), stats = S.stats; end
  if isfield(S, 'diagnostics'), diagnostics = S.diagnostics; end
  if isfield(S, 'info'), info = S.info; end

  [rMax, keepChains, rHatByParam] = maxRhatOverMonitoredParams( ...
    chains, monitorParams, p.Results.keepChainsMin, p.Results.rhatCritical);
  attemptInfo.rMax = rMax;
  attemptInfo.rHatByParam = rHatByParam;
  attemptInfo.keepChains = keepChains;
  fprintf('R-hat max=%.3f over {%s}; keep %d / %d chains\n', ...
    rMax, strjoin(fieldnames(rHatByParam), ', '), numel(keepChains), p.Results.nChains);

  if rMax <= p.Results.rhatCritical && numel(keepChains) >= p.Results.keepChainsMin
    fprintf('Stored chains meet R-hat <= %.3g — using saved fit\n', p.Results.rhatCritical);
    chains = subsetChainsHierarchical(chains, keepChains);
    showFinalHierarchicalFigures(chains, monitorParams, data, d, modelName, figuresDir, p.Results.saveFigures);
    converged = true;
    return;
  end
  fprintf('Stored chains do not meet R-hat <= %.3g — refitting\n', p.Results.rhatCritical);
end

if ~isfolder(storageDir)
  mkdir(storageDir);
end
if p.Results.resetThin
  clearMcmcState(storageDir, modelName);
end

if p.Results.singleAttempt
  nThin = p.Results.nThin;
  nBurnin = p.Results.nBurnin;
else
  [nThin, nBurnin] = initialMcmcFromState(storageDir, modelName, p.Results.nThin, p.Results.nBurnin);
end
converged = false;

while ~converged
  [stats, chains, diagnostics, info, rMax, keepChains, rHatByParam, nThin] = ...
    fitHierarchicalExecutionOnce(engine, modelName, monitorParams, data, initGenerator, ...
    nThin, nBurnin, p.Results, figuresDir);

  attemptInfo.rMax = rMax;
  attemptInfo.rHatByParam = rHatByParam;
  attemptInfo.nThin = nThin;
  attemptInfo.keepChains = keepChains;
  attemptInfo.chains = chains;
  try
    [thetaRMax, thetaRHatByP, thetaKeep] = assessHierarchicalThetaRhat( ...
      chains, data, d, modelName, p.Results.keepChainsMin, p.Results.rhatCritical);
    attemptInfo.thetaRMax = thetaRMax;
    attemptInfo.thetaRHatByParticipant = thetaRHatByP;
    attemptInfo.thetaKeepChains = thetaKeep;
  catch ME
    fprintf('theta R-hat assessment skipped: %s\n', ME.message);
  end

  if rMax <= p.Results.rhatCritical && numel(keepChains) >= p.Results.keepChainsMin
    chains = subsetChainsHierarchical(chains, keepChains);
    attemptInfo.chains = chains;
    converged = true;
    if p.Results.saveOnConverge
      save(storagePath, 'chains', 'stats', 'diagnostics', 'info', 'rHatByParam', '-v7.3');
      fprintf('Saved %s\n', storagePath);
      clearMcmcState(storageDir, modelName);
      showFinalHierarchicalFigures(chains, monitorParams, data, d, modelName, ...
        figuresDir, p.Results.saveFigures);
    else
      fprintf('Param R-hat ok but saveOnConverge=false — caller will decide\n');
    end
  else
    saveUnsuccessfulMcmcState(storageDir, modelName, nThin);
    if p.Results.singleAttempt
      fprintf('Not converged at thin=%d (max R-hat=%.3f)\n', nThin, rMax);
      break;
    end
    nThin = nextThin(nThin);
    fprintf('Not converged — next: thin=%d, burn-in=%g (unchanged)\n', nThin, nBurnin);
  end

  if converged || p.Results.singleAttempt
    if converged || p.Results.plotFailedRuns
      grtable(chains, p.Results.rhatCritical);
      codatable(chains);
    end
    break;
  end

  grtable(chains, p.Results.rhatCritical);
  codatable(chains);
end

end

function [stats, chains, diagnostics, info, rMax, keepChains, rHatByParam, nThin] = ...
    fitHierarchicalExecutionOnce(engine, modelName, monitorParams, data, initGenerator, ...
    nThin, nBurnin, opts, figuresDir)

tic;
modelPath = resolveJagsModelFile(fileparts(mfilename('fullpath')), ...
  sprintf('%s_%s.txt', modelName, engine));
[stats, chains, diagnostics, info] = callbayes(engine, ...
  'model', modelPath, ...
  'data', data, ...
  'outputname', 'samples', ...
  'init', initGenerator, ...
  'datafilename', modelName, ...
  'initfilename', modelName, ...
  'scriptfilename', modelName, ...
  'logfilename', fullfile('tmp', modelName), ...
  'nchains', opts.nChains, ...
  'nburnin', nBurnin, ...
  'nsamples', opts.nSamples, ...
  'monitorparams', monitorParams, ...
  'thin', nThin, ...
  'workingdir', fullfile('tmp', modelName), ...
  'verbosity', 0, ...
  'saveoutput', true, ...
  'allowunderscores', 1, ...
  'parallel', opts.doParallel);
fprintf('%s took %.1f s (thin=%d, burn-in=%g)\n', upper(engine), toc, nThin, nBurnin);

[rMax, keepChains, rHatByParam] = maxRhatOverMonitoredParams( ...
  chains, monitorParams, opts.keepChainsMin, opts.rhatCritical);
fprintf('R-hat max=%.3f over {%s}; keep %d / %d chains\n', ...
  rMax, strjoin(fieldnames(rHatByParam), ', '), numel(keepChains), opts.nChains);

if opts.plotFailedRuns
  plotParticipantParameterCIs(chains, monitorParams, ...
    'sgtitle', sprintf('%s | thin=%d | max R-hat=%.3f', modelName, nThin, rMax), ...
    'savePath', hierarchicalFigurePath(figuresDir, modelName, sprintf('thin%d', nThin), opts.saveFigures));
end

end

function showFinalHierarchicalFigures(chains, monitorParams, data, d, modelName, figuresDir, saveFigures)
plotParticipantParameterCIs(chains, monitorParams, ...
  'sgtitle', sprintf('%s | final', modelName), ...
  'savePath', hierarchicalFigurePath(figuresDir, modelName, 'final', saveFigures));
plotHierarchicalPosteriorPredictiveFigure(chains, data, d, modelName, ...
  'sgtitle', sprintf('%s | posterior predictive | final', modelName), ...
  'savePath', hierarchicalFigurePath(figuresDir, modelName, 'postPredFinal', saveFigures));
printParameterMeans(chains, monitorParams);
end

function printParameterMeans(chains, monitorParams)
for k = 1:numel(monitorParams)
  pname = monitorParams{k};
  mu = get_matrix_from_coda(chains, pname, @mean);
  fprintf('%s mean (p1..p5): ', pname);
  fprintf('%.4g ', mu(1:min(5, numel(mu))));
  fprintf('\n');
end
end

function chainsOut = subsetChainsHierarchical(chainsIn, keepChains)
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
