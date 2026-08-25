function runHierarchicalModelStan(modelName, monitorParams, initGenerator, varargin)
%RUNHIERARCHICALMODELSTAN  Hierarchical model via Stan (Trinity / CmdStan).
%
%   Mirrors runHierarchicalExecutionModel.m but uses engine 'stan' and
%   models/stan/{modelName}.stan. Chains and figures go to stan/storage and
%   stan/figures. Reuses parent-folder diagnostics and plotting utilities.

p = inputParser;
addParameter(p, 'preLoad', true, @islogical);
addParameter(p, 'rhatCritical', 1.1, @isnumeric);
addParameter(p, 'keepChainsMin', 4, @isnumeric);
addParameter(p, 'nChains', 8, @isnumeric);
addParameter(p, 'nBurnin', 1e3, @isnumeric);
addParameter(p, 'nSamples', 1e3, @isnumeric);
addParameter(p, 'nThin', 1, @isnumeric);
addParameter(p, 'doParallel', true, @islogical);
addParameter(p, 'dataName', 'intertemporalChoice', @ischar);
addParameter(p, 'dataDir', '', @ischar);
addParameter(p, 'saveFigures', true, @islogical);
addParameter(p, 'resetThin', false, @islogical);
parse(p, varargin{:});

engine = 'stan';
monitorParams = cellstr(monitorParams(:));
stanDir = fileparts(mfilename('fullpath'));
modelsDir = fileparts(stanDir);
generalDir = fullfile(modelsDir, '..', 'general');
figuresDir = fullfile(stanDir, 'figures');
storageDir = fullfile(stanDir, 'storage');
sourceModelPath = fullfile(stanDir, [modelName '.stan']);
workingDir = fullfile(stanDir, 'tmp', modelName);

addpath(modelsDir);
addpath(stanDir);
addpath(generalDir);
cleanupObj = onCleanup(@() rmpath(generalDir));

[data, d] = prepareIntertemporalChoiceData(p.Results.dataName, p.Results.dataDir);

if nargin < 3 || isempty(initGenerator)
  if any(strcmp(monitorParams, 'sigma'))
    initGenerator = @() struct('sigma', rand(d.nParticipants, 1) * 0.4 + 0.05);
  else
    initGenerator = @() struct('epsilon', rand(d.nParticipants, 1) * 0.4 + 0.05);
  end
end

if ~isfile(sourceModelPath)
  error('runHierarchicalModelStan:missingModel', 'Stan model not found: %s', sourceModelPath);
end
modelPath = stageStanModelForCmdStan(stanDir, modelName, workingDir);

fileName = sprintf('%s_%s_%s.mat', modelName, p.Results.dataName, engine);
storagePath = fullfile(storageDir, fileName);

if p.Results.preLoad && isfile(storagePath)
  fprintf('Loading stored chains: %s\n', storagePath);
  S = load(storagePath, 'chains', 'stats', 'diagnostics', 'info');
  chains = S.chains;
  if isfield(S, 'stats'), stats = S.stats; end
  if isfield(S, 'diagnostics'), diagnostics = S.diagnostics; end
  if isfield(S, 'info'), info = S.info; end
  showFinalHierarchicalFiguresStan(chains, monitorParams, data, d, modelName, figuresDir, p.Results.saveFigures);
  return;
end

if ~isfolder(storageDir)
  mkdir(storageDir);
end
if p.Results.resetThin
  clearMcmcState(storageDir, modelName);
end
[nThin, nBurnin] = initialMcmcFromState(storageDir, modelName, p.Results.nThin, p.Results.nBurnin);
converged = false;

while ~converged
  tic;
  [stats, chains, diagnostics, info] = callbayes(engine, ...
    'model', modelPath, ...
    'data', data, ...
    'outputname', 'samples', ...
    'init', initGenerator, ...
    'datafilename', modelName, ...
    'initfilename', modelName, ...
    'scriptfilename', modelName, ...
    'logfilename', fullfile('tmp', modelName), ...
    'nchains', p.Results.nChains, ...
    'nburnin', nBurnin, ...
    'nsamples', p.Results.nSamples, ...
    'monitorparams', monitorParams, ...
    'thin', nThin, ...
    'workingdir', workingDir, ...
    'verbosity', 0, ...
    'saveoutput', true, ...
    'parallel', p.Results.doParallel, ...
    'allowunderscores', true);
  fprintf('%s took %.1f s (thin=%d, burn-in=%g)\n', upper(engine), toc, nThin, nBurnin);

  [rMax, keepChains, rHatByParam] = maxRhatOverMonitoredParams( ...
    chains, monitorParams, p.Results.keepChainsMin, p.Results.rhatCritical);
  fprintf('R-hat max=%.3f over {%s}; keep %d / %d chains\n', ...
    rMax, strjoin(fieldnames(rHatByParam), ', '), numel(keepChains), p.Results.nChains);

  plotParticipantParameterCIs(chains, monitorParams, ...
    'sgtitle', sprintf('%s | thin=%d | max R-hat=%.3f', modelName, nThin, rMax), ...
    'savePath', hierarchicalFigurePath(figuresDir, modelName, sprintf('thin%d', nThin), p.Results.saveFigures));

  if rMax < p.Results.rhatCritical && numel(keepChains) >= p.Results.keepChainsMin
    chains = subsetChainsHierarchical(chains, keepChains);
    converged = true;
    save(storagePath, 'chains', 'stats', 'diagnostics', 'info', 'rHatByParam', '-v7.3');
    fprintf('Saved %s\n', storagePath);
    clearMcmcState(storageDir, modelName);
  else
    saveUnsuccessfulMcmcState(storageDir, modelName, nThin, nBurnin);
    nThinNext = nextThin(nThin);
    nBurninNext = nextBurnin(nBurnin);
    fprintf('Not converged — next: thin=%d (was %d), burn-in=%g (was %g)\n', ...
      nThinNext, nThin, nBurninNext, nBurnin);
    nThin = nThinNext;
    nBurnin = nBurninNext;
  end

  grtable(chains, p.Results.rhatCritical);
  codatable(chains);
end

showFinalHierarchicalFiguresStan(chains, monitorParams, data, d, modelName, figuresDir, p.Results.saveFigures);

end

function showFinalHierarchicalFiguresStan(chains, monitorParams, data, d, modelName, figuresDir, saveFigures)
plotParticipantParameterCIs(chains, monitorParams, ...
  'sgtitle', sprintf('%s | final', modelName), ...
  'savePath', hierarchicalFigurePath(figuresDir, modelName, 'final', saveFigures));
plotHierarchicalPosteriorPredictiveFigure(chains, data, d, modelName, ...
  'sgtitle', sprintf('%s | posterior predictive | final', modelName), ...
  'savePath', hierarchicalFigurePath(figuresDir, modelName, 'postPredFinal', saveFigures));
printParameterMeansStan(chains, monitorParams);
end

function printParameterMeansStan(chains, monitorParams)
for k = 1:numel(monitorParams)
  pname = monitorParams{k};
  mu = get_matrix_from_coda(chains, pname, @mean);
  fprintf('%s mean (p1..p5): ', pname);
  fprintf('%.4g ', mu(1:min(5, numel(mu))));
  fprintf('\n');
end
end

function modelPath = stageStanModelForCmdStan(stanDir, modelName, workingDir)
%STAGESTANMODELFORCMDSTAN  Merge common_functions into model for CmdStan (no #include).
if ~isfolder(workingDir)
  mkdir(workingDir);
end
commonPath = fullfile(stanDir, 'common_functions.stan');
modelPath = fullfile(workingDir, [modelName '.stan']);
commonText = fileread(commonPath);
modelText = fileread(fullfile(stanDir, [modelName '.stan']));
modelText = regexprep(modelText, '#include\s+"common_functions\.stan"\s*', '', 'once');
fid = fopen(modelPath, 'w');
assert(fid > 0, 'Cannot write staged Stan model: %s', modelPath);
fprintf(fid, '%s\n%s', commonText, modelText);
fclose(fid);
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
