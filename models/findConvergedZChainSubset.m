function [ok, keepChains, rMax, rHatEach] = findConvergedZChainSubset( ...
  zAcc, subsetSizeMin, rhatCritical)
%FINDCONVERGEDZCHAINSUBSET  Largest z-chain subset with all R-hat <= threshold.
%
%   Searches combinations of size nChains, nChains-1, ..., subsetSizeMin.
%   Returns the largest size that has at least one combo with max participant
%   z R-hat <= rhatCritical (among ties at that size, the lowest max R-hat).
%
%   If none qualify, ok=false and keepChains/rHatEach refer to all chains.
%
%   See also: MAXRHATOVERPARTICIPANTS, GELMANRUBIN

if nargin < 3 || isempty(rhatCritical)
  rhatCritical = 1.1;
end
if nargin < 2 || isempty(subsetSizeMin)
  subsetSizeMin = 4;
end

mats = codaIndexedMatrices(zAcc, 'z');
nP = numel(mats);
if nP < 1
  ok = false;
  keepChains = [];
  rMax = nan;
  rHatEach = [];
  return;
end
nC = size(mats{1}, 2);
subsetSizeMin = max(2, min(subsetSizeMin, nC));

for k = nC:-1:subsetSizeMin
  combos = nchoosek(1:nC, k);
  bestMax = inf;
  bestCombo = [];
  bestEach = [];
  for i = 1:size(combos, 1)
    ch = combos(i, :);
    rEach = nan(nP, 1);
    for p = 1:nP
      rEach(p) = gelmanRubinSafeLocal(mats{p}(:, ch));
    end
    thisMax = max(rEach);
    if thisMax < bestMax
      bestMax = thisMax;
      bestCombo = ch;
      bestEach = rEach;
    end
  end
  if bestMax <= rhatCritical
    ok = true;
    keepChains = bestCombo;
    rMax = bestMax;
    rHatEach = bestEach;
    return;
  end
end

% No qualifying subset: report all-chain diagnostics.
ok = false;
keepChains = 1:nC;
rHatEach = nan(nP, 1);
for p = 1:nP
  rHatEach(p) = gelmanRubinSafeLocal(mats{p});
end
rMax = max(rHatEach);
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
