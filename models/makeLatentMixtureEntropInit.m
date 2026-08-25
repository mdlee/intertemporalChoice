function initGenerator = makeLatentMixtureEntropInit(storageDir, dataName, engine, nParticipants, varargin)
%MAKELATENTMIXTUREENTROPINIT  Warm-start latent-mixture entrop from hierarchical fits.
%
%   Loads *ExecutionHierarchical_entrop (and SS/LL) mats under storageDir,
%   maps fields to mixture participant names, sets Trinity-safe group
%   hyperpriors muFoo / tauFoo from across-participant means / precisions,
%   and returns an init function handle for callbayes.
%
%   Name-value:
%     zMode       ([])     — nParticipants x 1 categorical z; random if empty
%     jitterFrac  (0.05)   — relative Gaussian jitter per chain
%     seed        ([])     — optional RNG seed
%     nModels     (11)     — for random z draws
%     summaryFn   (@median)— coda summary for participant / alpha inits
%                            (group mu*/tau* derived from those summaries)

p = inputParser;
addParameter(p, 'zMode', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'jitterFrac', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
addParameter(p, 'nModels', 11, @isnumeric);
addParameter(p, 'summaryFn', @median, @(f) isa(f, 'function_handle'));
parse(p, varargin{:});

nParticipants = double(nParticipants);
jitterFrac = p.Results.jitterFrac;
nModels = double(p.Results.nModels);
summaryFn = p.Results.summaryFn;

% Each row: {hierarchicalStem, { {hierField, mixtureField}, ... }}
cognitiveSpecs = { ...
  'exponentialExecutionHierarchical_entrop', ...
    {{'kappa', 'kappaEX'}}; ...
  'hyperbolicExecutionHierarchical_entrop', ...
    {{'kappa', 'kappaHC'}}; ...
  'hyperboloidExecutionHierarchical_entrop', ...
    {{'kappa', 'kappaHY'}, {'tau', 'tauHY'}}; ...
  'proportionalDifferencesExecutionHierarchical_entrop', ...
    {{'delta', 'deltaPD'}}; ...
  'directDifferencesExecutionHierarchical_entrop', ...
    {{'delta', 'deltaDD'}, {'omega', 'omegaDD'}}; ...
  'tradeoffExecutionHierarchical_entrop', ...
    {{'gamma', 'gammaTR'}, {'tau', 'tauTR'}, {'kappa', 'kappaTR'}, {'vartheta', 'varthetaTR'}}; ...
  'unifiedTradeoffExecutionHierarchical_entrop', ...
    {{'gamma', 'gammaUT'}, {'tau', 'tauUT'}, {'kappa', 'kappaUT'}, ...
     {'vartheta', 'varthetaUT'}, {'eta', 'etaUT'}}; ...
  'itchExecutionHierarchical_entrop', ...
    {{'beta0', 'beta0IT'}, {'betaRA', 'betaRAIT'}, {'betaRR', 'betaRRIT'}, ...
     {'betaTA', 'betaTAIT'}, {'betaTR', 'betaTRIT'}} ...
  };

means = struct();
wSum = zeros(nParticipants, 1);
wCount = 0;

for s = 1:size(cognitiveSpecs, 1)
  stem = cognitiveSpecs{s, 1};
  fieldMap = cognitiveSpecs{s, 2};
  path = storageMatPath(storageDir, stem, dataName, engine, '');
  if ~isfile(path)
    error('makeLatentMixtureEntropInit:missingFit', ...
      'Missing hierarchical fit for warm-start: %s', path);
  end
  fprintf('Init source (%s): %s\n', func2str(summaryFn), path);
  S = load(path, 'chains');
  for f = 1:numel(fieldMap)
    hierName = fieldMap{f}{1};
    mixName = fieldMap{f}{2};
    mu = get_matrix_from_coda(S.chains, hierName, summaryFn);
    mu = mu(:);
    if numel(mu) < nParticipants
      error('makeLatentMixtureEntropInit:badChains', ...
        '%s in %s has %d elements; need %d.', hierName, stem, numel(mu), nParticipants);
    end
    means.(mixName) = mu(1:nParticipants);
  end
  if codaHasParam(S.chains, 'w')
    ww = get_matrix_from_coda(S.chains, 'w', summaryFn);
    ww = ww(:);
    wSum = wSum + ww(1:nParticipants);
    wCount = wCount + 1;
  end
end

if wCount == 0
  means.w = ones(nParticipants, 1) * 2;
else
  means.w = wSum / wCount;
end

means.alphaLL = loadAlphaVector(storageDir, 'LL', dataName, engine, nParticipants, summaryFn);
means.alphaSS = loadAlphaVector(storageDir, 'SS', dataName, engine, nParticipants, summaryFn);

if isempty(p.Results.zMode)
  z0 = [];
else
  z0 = round(p.Results.zMode(:));
  z0 = max(1, min(nModels, z0(1:nParticipants)));
end

% Group hyperpriors from participant inits (Trinity-safe camelCase names)
groupMeans = struct();
groupTaus = struct();
partNames = fieldnames(means);
for k = 1:numel(partNames)
  pname = partNames{k};
  if any(strcmp(pname, {'alphaLL', 'alphaSS'}))
    continue;
  end
  v = means.(pname);
  muName = ['mu' upper(pname(1)) pname(2:end)];
  tauName = ['tau' upper(pname(1)) pname(2:end)];
  % Special-case single-letter w → muW / tauW
  if strcmp(pname, 'w')
    muName = 'muW';
    tauName = 'tauW';
  end
  groupMeans.(muName) = median(v);
  s2 = var(v, 0);
  if ~(isfinite(s2) && s2 > 1e-8)
    s2 = 0.25;
  end
  prec = 1 / s2;
  prec = min(max(prec, 0.05 + 1e-6), 25 - 1e-6);
  groupTaus.(tauName) = prec;
end

if ~isempty(p.Results.seed)
  rng(p.Results.seed);
end

initGenerator = @() drawMixtureInit(means, groupMeans, groupTaus, z0, ...
  nParticipants, nModels, jitterFrac);
end

function s = drawMixtureInit(means, groupMeans, groupTaus, z0, nP, nModels, jitterFrac)
s = struct();
if isempty(z0)
  s.z = randi(nModels, nP, 1);
else
  s.z = z0;
end

partNames = fieldnames(means);
for k = 1:numel(partNames)
  pname = partNames{k};
  base = means.(pname);
  [lo, hi] = mixtureParamBounds(pname);
  jittered = base + jitterFrac .* max(abs(base), 0.05) .* randn(size(base));
  s.(pname) = clampVec(jittered, lo, hi);
end

muNames = fieldnames(groupMeans);
for k = 1:numel(muNames)
  mname = muNames{k};
  base = groupMeans.(mname);
  % Infer participant name for bounds: muKappaEX → kappaEX, muW → w
  if strcmp(mname, 'muW')
    pname = 'w';
  else
    pname = [lower(mname(3)) mname(4:end)];
  end
  [lo, hi] = mixtureParamBounds(pname);
  jittered = base + jitterFrac .* max(abs(base), 0.05) .* randn();
  s.(mname) = clampVec(jittered, lo, hi);
end

tauNames = fieldnames(groupTaus);
for k = 1:numel(tauNames)
  tname = tauNames{k};
  base = groupTaus.(tname);
  jittered = base * exp(jitterFrac * randn());
  s.(tname) = clampVec(jittered, 0.05, 25);
end
end

function alpha = loadAlphaVector(storageDir, modelName, dataName, engine, nP, summaryFn)
if nargin < 6 || isempty(summaryFn)
  summaryFn = @median;
end
path = storageMatPath(storageDir, modelName, dataName, engine, '');
if ~isfile(path)
  warning('makeLatentMixtureEntropInit:missingAlpha', ...
    'Missing %s; using dbeta-ish default alpha.', path);
  alpha = 0.95 * ones(nP, 1);
  return;
end
S = load(path, 'chains');
alpha = get_matrix_from_coda(S.chains, 'alpha', summaryFn);
alpha = alpha(:);
if numel(alpha) < nP
  error('makeLatentMixtureEntropInit:badAlpha', ...
    '%s alpha has %d elements; need %d.', modelName, numel(alpha), nP);
end
alpha = alpha(1:nP);
end

function tf = codaHasParam(chains, pname)
fn = fieldnames(chains);
tf = any(strcmp(fn, pname)) || any(strncmp(fn, [pname '_'], numel(pname) + 1));
end

function [lo, hi] = mixtureParamBounds(pname)
switch pname
  case {'kappaEX', 'kappaHC', 'kappaHY', 'tauHY', ...
      'gammaTR', 'tauTR', 'kappaTR', 'varthetaTR', ...
      'gammaUT', 'tauUT', 'kappaUT', 'varthetaUT'}
    lo = 1e-4; hi = 100;
  case 'etaUT'
    lo = 0; hi = 500;
  case 'deltaPD'
    lo = -2; hi = 2;
  case 'deltaDD'
    lo = -2000; hi = 2000;
  case 'omegaDD'
    lo = 1e-6; hi = 1 - 1e-6;
  case 'beta0IT'
    lo = -2; hi = 2;
  case {'betaRAIT', 'betaRRIT'}
    lo = 0; hi = 2;
  case {'betaTAIT', 'betaTRIT'}
    lo = -2; hi = 0;
  case 'w'
    lo = 0; hi = 20;
  case {'alphaLL', 'alphaSS'}
    lo = 1e-6; hi = 1 - 1e-6;
  otherwise
    lo = -inf; hi = inf;
end
end

function x = clampVec(x, lo, hi)
x = min(max(x, lo), hi);
epsIn = 1e-6;
if isfinite(lo)
  x = max(x, lo + epsIn);
end
if isfinite(hi)
  x = min(x, hi - epsIn);
end
x = x(:);
end
