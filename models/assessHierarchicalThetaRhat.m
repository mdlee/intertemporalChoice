function [rMax, rHatByParticipant, keepChains] = assessHierarchicalThetaRhat( ...
  chains, data, d, modelName, keepChainsMin, rhatCritical)
%ASSESSHIERARCHICALTHETARHAT  Gelman-Rubin R-hat on reconstructed theta (PP probs).
%
%   JAGS defines theta[i,j] but does not monitor it (no PP replicates either).
%   This rebuilds trial-mean P(LL) per MCMC draw — the quantity that drives the
%   MATLAB posterior-predictive figures — and reports R-hat by participant.

if nargin < 5 || isempty(keepChainsMin)
  keepChainsMin = 8;
end
if nargin < 6 || isempty(rhatCritical)
  rhatCritical = 1.4;
end

nP = double(d.nParticipants);
rHatByParticipant = nan(nP, 1);
keepChains = [];

for pp = 1:nP
  ThMean = trialMeanThetaMatrix(chains, data, modelName, pp);  % nSamples × nChains
  rHatByParticipant(pp) = gelmanRubinSafeMatrix(ThMean);

  if rHatByParticipant(pp) <= rhatCritical
    kc = 1:size(ThMean, 2);
  else
    [kc, ~] = findKeepChains(ThMean, keepChainsMin, rhatCritical);
    if isempty(kc)
      kc = 1:size(ThMean, 2);
    end
  end

  if pp == 1
    keepChains = kc;
  else
    keepChains = intersect(keepChains, kc, 'stable');
  end
end

rMax = max(rHatByParticipant);
if isempty(keepChains)
  keepChains = 1:size(trialMeanThetaMatrix(chains, data, modelName, 1), 2);
end

fprintf('theta (trial-mean P(LL) / PP) R-hat: max=%.3f, min=%.3f; keep %d chains\n', ...
  rMax, min(rHatByParticipant), numel(keepChains));
[sorted, idx] = sort(rHatByParticipant, 'descend');
for k = 1:min(5, numel(sorted))
  fprintf('  p%d theta R-hat=%.3f\n', idx(k), sorted(k));
end
end

function ThMean = trialMeanThetaMatrix(chains, data, modelName, pp)
paramNames = structuralParamsForModel(modelName);
nParams = numel(paramNames);
mats = cell(nParams, 1);
for k = 1:nParams
  allP = codaIndexedMatrices(chains, paramNames{k});
  mats{k} = allP{pp};  % nSamples × nChains
end
nS = size(mats{1}, 1);
nC = size(mats{1}, 2);

rLL = data.rLL(pp, :);
rSS = data.rSS(pp, :);
tLL = data.tLL(pp, :);
tSS = data.tSS(pp, :);

ThMean = nan(nS, nC);
for c = 1:nC
  cols = cell(nParams, 1);
  for k = 1:nParams
    cols{k} = mats{k}(:, c);
  end
  Th = thetaFromColumns(cols, paramNames, modelName, rLL, rSS, tLL, tSS);
  ThMean(:, c) = mean(Th, 2);
end
end

function Th = thetaFromColumns(cols, paramNames, modelName, rLL, rSS, tLL, tSS)
useEntrop = endsWith(modelName, '_entrop');
baseName = modelName;
if useEntrop
  baseName = extractBefore(modelName, '_entrop');
end
if ~strcmp(baseName, 'unifiedTradeoffExecutionHierarchical')
  error('assessHierarchicalThetaRhat:unsupported', 'Unhandled model %s', modelName);
end
% Expected order: gamma, tau, kappa, vartheta, eta, w|epsilon
gamma = cols{1}; tau = cols{2}; kappa = cols{3}; vartheta = cols{4};
eta = cols{5}; w = cols{6};
if ~useEntrop
  error('assessHierarchicalThetaRhat:unsupported', 'Non-entrop unused here');
end
Th = entropTradeoffThetaLocal(gamma, tau, kappa, vartheta, w, rLL, rSS, tLL, tSS, true, eta);
Th = min(max(Th, 0), 1);
end

function names = structuralParamsForModel(modelName)
useEntrop = endsWith(modelName, '_entrop');
if useEntrop
  names = {'gamma', 'tau', 'kappa', 'vartheta', 'eta', 'w'};
else
  names = {'gamma', 'tau', 'kappa', 'vartheta', 'eta', 'epsilon'};
end
end

function Th = entropTradeoffThetaLocal(gamma, tau, kappa, vartheta, w, rLL, rSS, tLL, tSS, unified, eta)
nD = numel(gamma);
nT = numel(rLL);
Th = nan(nD, nT);
for s = 1:nD
  g = gamma(s); tv = tau(s); k = kappa(s); vt = vartheta(s);
  if unified
    vLL = (1 / g) * log(1 + g * max(rLL - eta(s), 1e-9));
  else
    vLL = (1 / g) * log(1 + g * rLL);
  end
  vSS = (1 / g) * log(1 + g * rSS);
  wLL = (1 / tv) * log(1 + tv * tLL);
  wSS = (1 / tv) * log(1 + tv * tSS);
  qV = vLL - vSS;
  dw = wLL - wSS;
  qT = k * log(1 + max((dw / vt) .^ vt, 0));
  pick = (qV - qT) > 0;
  pLL = 1 / (1 + exp(-w(s)));
  pSS = 1 / (1 + exp(w(s)));
  Th(s, pick) = pLL;
  Th(s, ~pick) = pSS;
end
end

function rhat = gelmanRubinSafeMatrix(x)
x = double(x);
if isempty(x)
  rhat = 1;
  return;
end
if isscalar(unique(round(x(:), 8)))
  rhat = 1;
  return;
end
nChains = size(x, 2);
if nChains < 2
  rhat = 1;
  return;
end
w = sum(var(x, 0, 1)) / nChains;
if w <= 0 || ~isfinite(w)
  rhat = 1;
  return;
end
rhat = gelmanrubin(x, 0, 1, 'rhat');
if ~isfinite(rhat)
  rhat = 1;
end
end
