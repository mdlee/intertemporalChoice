function [ok, keepChains, rMax, rMaxAll] = findConvergedMonitoredChainSubset( ...
  coda, paramNames, subsetSizeMin, rhatCritical, participantIndices)
%FINDCONVERGEDMONITOREDCHAINSUBSET  Largest chain subset with all R-hat <= gate.
%
%   Searches combinations of size nChains, nChains-1, ..., subsetSizeMin.
%   A subset passes if every monitored parameter (every selected participant
%   index) has Gelman-Rubin R-hat <= rhatCritical.
%   Among passing subsets, the largest size wins; ties take the lowest max R-hat.
%   If none pass, keepChains is the best subset of size subsetSizeMin (lowest
%   max R-hat) so callers can track "best 4-chain R-hat" even while failing.
%
%   See also: FINDCONVERGEDZCHAINSUBSET, MAXRHATOVERMONITOREDPARAMS

if nargin < 5
  participantIndices = [];
end
if nargin < 4 || isempty(rhatCritical)
  rhatCritical = 1.05;
end
if nargin < 3 || isempty(subsetSizeMin)
  subsetSizeMin = 4;
end

paramNames = cellstr(paramNames(:));
matsByParam = cell(numel(paramNames), 1);
nC = [];
for p = 1:numel(paramNames)
  mats = codaIndexedMatrices(coda, paramNames{p});
  matsByParam{p} = mats;
  if isempty(nC)
    nC = size(mats{1}, 2);
  end
end
if isempty(nC) || nC < 2
  ok = false;
  keepChains = 1:max(1, double(nC));
  rMax = nan;
  rMaxAll = nan;
  return;
end

subsetSizeMin = max(2, min(subsetSizeMin, nC));
rMaxAll = comboMaxRhat(matsByParam, 1:nC, participantIndices);

bestFailMax = inf;
bestFailCombo = 1:nC;
for k = nC:-1:subsetSizeMin
  combos = nchoosek(1:nC, k);
  bestMax = inf;
  bestCombo = [];
  for i = 1:size(combos, 1)
    ch = combos(i, :);
    thisMax = comboMaxRhat(matsByParam, ch, participantIndices);
    if thisMax < bestMax
      bestMax = thisMax;
      bestCombo = ch;
    end
  end
  if bestMax <= rhatCritical
    ok = true;
    keepChains = bestCombo;
    rMax = bestMax;
    return;
  end
  if k == subsetSizeMin && bestMax < bestFailMax
    bestFailMax = bestMax;
    bestFailCombo = bestCombo;
  end
end

ok = false;
keepChains = bestFailCombo;
rMax = bestFailMax;
end

function rMax = comboMaxRhat(matsByParam, ch, participantIndices)
rMax = -inf;
for p = 1:numel(matsByParam)
  mats = matsByParam{p};
  nElem = numel(mats);
  if isempty(participantIndices)
    idx = 1:nElem;
  else
    idx = unique(double(participantIndices(:)'), 'stable');
    idx = idx(idx >= 1 & idx <= nElem);
  end
  for j = idx
    rhat = gelmanRubinSafeLocal(mats{j}(:, ch));
    rMax = max(rMax, rhat);
  end
end
if rMax < 0
  rMax = nan;
end
end

function rhat = gelmanRubinSafeLocal(x)
x = double(x);
if isempty(x)
  rhat = 1;
  return;
end
if isscalar(unique(round(x(:))))
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
