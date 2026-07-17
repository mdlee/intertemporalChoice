function [rMax, keepChains, rHatEach] = maxRhatOverParticipants(coda, minChains, desiredRhat, paramName)
%MAXRHATOVERPARTICIPANTS  Worst Gelman-Rubin R-hat over an indexed CODA parameter.
%
%   RMAX = MAXRHATOVERPARTICIPANTS(CODA) — worst R-hat across z_1, z_2, ...
%
%   [RMAX, KEEPCHAINS, RHATEACH] = MAXRHATOVERPARTICIPANTS(CODA, MINCHAINS, DESIRED)
%   RHATEACH is a vector of per-index R-hats. If all chains share the same value
%   (common for discrete z), R-hat is treated as 1 instead of Inf.
%
%   See also: CODAINDEXEDMATRICES, FINDKEEPCHAINS, GELMANRUBIN

if nargin < 4 || isempty(paramName)
  paramName = 'z';
end
if nargin < 3 || isempty(desiredRhat)
  desiredRhat = inf;
end
if nargin < 2 || isempty(minChains)
  minChains = 2;
end

matrices = codaIndexedMatrices(coda, paramName);
nElem = numel(matrices);
rHatEach = nan(nElem, 1);
rMax = -inf;
keepChains = [];

for k = 1:nElem
  x = matrices{k};
  rHatEach(k) = gelmanRubinSafe(x);
  rMax = max(rMax, rHatEach(k));

  if rHatEach(k) <= desiredRhat
    kc = 1:size(x, 2);
  else
    [kc, ~] = findKeepChains(x, minChains, desiredRhat);
    if isempty(kc)
      kc = 1:size(x, 2);
    end
  end

  if k == 1
    keepChains = kc;
  else
    keepChains = intersect(keepChains, kc, 'stable');
  end
end

if isempty(keepChains) && nElem > 0
  keepChains = 1:size(matrices{1}, 2);
end

end

function rhat = gelmanRubinSafe(x)
% Constant chains (e.g. stuck discrete z) => Rhat 1; non-finite => 1.

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
