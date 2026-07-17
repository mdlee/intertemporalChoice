function [rMax, keepChains, rHatByParam] = maxRhatOverMonitoredParams(chains, monitorParams, minChains, desiredRhat)
%MAXRHATOVERMONITOREDPARAMS  Worst R-hat across a cell list of indexed parameters.

if nargin < 4 || isempty(desiredRhat)
  desiredRhat = inf;
end
if nargin < 3 || isempty(minChains)
  minChains = 2;
end

monitorParams = cellstr(monitorParams(:));
nPar = numel(monitorParams);
rHatByParam = struct();
rMax = -inf;
keepChains = [];

for p = 1:nPar
  pname = monitorParams{p};
  [rHat, kc] = maxRhatOverParticipants(chains, minChains, desiredRhat, pname);
  rHatByParam.(pname) = rHat;
  rMax = max(rMax, rHat);
  if p == 1
    keepChains = kc;
  else
    keepChains = intersect(keepChains, kc, 'stable');
  end
end

if isempty(keepChains)
  mats = codaIndexedMatrices(chains, monitorParams{1});
  if ~isempty(mats)
    keepChains = 1:size(mats{1}, 2);
  end
end

end
