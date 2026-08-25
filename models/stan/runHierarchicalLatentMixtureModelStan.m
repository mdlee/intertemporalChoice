function runHierarchicalLatentMixtureModelStan(modelName, initGenerator, varargin)
%RUNHIERARCHICALLATENTMIXTUREMODELSTAN  Latent mixture via Stan (Trinity / CmdStan).
%
%   Mirrors runHierarchicalLatentMixtureModel.m. Discrete z is marginalized in
%   Stan; generated quantities provide argmax-model z for plotting (see .stan).

p = inputParser;
addParameter(p, 'preLoad', false, @islogical);
addParameter(p, 'resetThin', false, @islogical);
addParameter(p, 'rhatCritical', 1.1, @isnumeric);
addParameter(p, 'keepChainsMin', 12, @isnumeric);
addParameter(p, 'nChains', 12, @isnumeric);
addParameter(p, 'nBurnin', 1e3, @isnumeric);
addParameter(p, 'nSamples', 2e3, @isnumeric);
addParameter(p, 'nThin', 1, @isnumeric);
addParameter(p, 'nModels', 11, @isnumeric);
addParameter(p, 'doParallel', true, @islogical);
addParameter(p, 'dataName', 'intertemporalChoice', @ischar);
addParameter(p, 'dataDir', '', @ischar);
addParameter(p, 'saveFigures', true, @islogical);
parse(p, varargin{:});

engine = 'stan';
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
data.nModels = p.Results.nModels;

if ~isfile(sourceModelPath)
  error('runHierarchicalLatentMixtureModelStan:missingModel', 'Stan model not found: %s', sourceModelPath);
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
  showFinalLatentMixtureFiguresStan(chains, data, d, modelName, figuresDir, p.Results.saveFigures);
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
    'monitorparams', latentMixtureHierarchicalMonitorParams(), ...
    'thin', nThin, ...
    'workingdir', workingDir, ...
    'verbosity', 0, ...
    'saveoutput', true, ...
    'parallel', p.Results.doParallel, ...
    'allowunderscores', true);
  fprintf('%s took %.1f s (thin=%d, burn-in=%g)\n', upper(engine), toc, nThin, nBurnin);

  [zRhat, keepChains, rHatEach] = maxRhatOverParticipants( ...
    chains, p.Results.keepChainsMin, p.Results.rhatCritical);
  fprintf('z R-hat: max=%.3f, min=%.3f, %d/%d at 1 (constant z); keep %d / %d chains\n', ...
    zRhat, min(rHatEach), sum(rHatEach <= 1 + 1e-9), numel(rHatEach), ...
    numel(keepChains), p.Results.nChains);
  plotLatentMixtureZPosterior(chains, ...
    'nModels', data.nModels, ...
    'sgtitle', sprintf('%s | thin=%d | max R-hat=%.3f', modelName, nThin, zRhat), ...
    'savePath', hierarchicalFigurePath(figuresDir, modelName, sprintf('thin%d', nThin), p.Results.saveFigures));

  if zRhat < p.Results.rhatCritical && numel(keepChains) >= p.Results.keepChainsMin
    chains = subsetChainsLatentMixtureStan(chains, keepChains);
    converged = true;
    save(storagePath, 'chains', 'stats', 'diagnostics', 'info', '-v7.3');
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

showFinalLatentMixtureFiguresStan(chains, data, d, modelName, figuresDir, p.Results.saveFigures);

end

function showFinalLatentMixtureFiguresStan(chains, data, d, modelName, figuresDir, saveFigures)
zMean = get_matrix_from_coda(chains, 'z', @mean);
zMode = get_matrix_from_coda(chains, 'z', @mode);
disp(table((1:d.nParticipants)', zMean(:), zMode(:), 'VariableNames', {'p', 'zMean', 'zMode'}));
plotLatentMixtureZPosterior(chains, ...
  'nModels', data.nModels, ...
  'sgtitle', sprintf('%s | final', modelName), ...
  'savePath', hierarchicalFigurePath(figuresDir, modelName, 'final', saveFigures));
end

function modelPath = stageStanModelForCmdStan(stanDir, modelName, workingDir)
if ~isfolder(workingDir)
  mkdir(workingDir);
end
commonText = fileread(fullfile(stanDir, 'common_functions.stan'));
modelPath = fullfile(workingDir, [modelName '.stan']);
modelText = fileread(fullfile(stanDir, [modelName '.stan']));
modelText = regexprep(modelText, '#include\s+"common_functions\.stan"\s*', '', 'once');
fid = fopen(modelPath, 'w');
assert(fid > 0, 'Cannot write staged Stan model: %s', modelPath);
fprintf(fid, '%s\n%s', commonText, modelText);
fclose(fid);
end

function chainsOut = subsetChainsLatentMixtureStan(chainsIn, keepChains)
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
