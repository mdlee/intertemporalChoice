function [logW, khat] = psisSmoothLogWeights(logRatios)
%PSISSMOOTHLOGWEIGHTS  Pareto-smoothed importance sampling (Vehtari et al.).
%
%   LOGRATIOS are log importance ratios, e.g. -log p(y_i | θ_s) for LOO.
%   Returns log smoothed weights (unnormalized) and the GPD shape khat.
%
%   Follows the PSIS algorithm used by Yao, Vehtari & Gelman (2022) / loo.

logRatios = logRatios(:);
S = numel(logRatios);
khat = 0;
if S < 5
  logW = logRatios;
  return;
end

[lrSorted, ord] = sort(logRatios);
logW = logRatios;

M = min(floor(S / 5), max(1, floor(3 * sqrt(S))));
cutoff = lrSorted(S - M);
tail = lrSorted(S - M + 1:S);
exc = tail - cutoff;
exc = max(exc, realmin);

% Method-of-moments GPD fit on excesses (Zhang & Stephens style MOM)
m = mean(exc);
v = var(exc, 1);
if ~(isfinite(m) && m > 0 && isfinite(v) && v > 0)
  return;
end
khat = (m * m / v - 1) / 2;
sigma = m * (1 - khat);
if ~(isfinite(khat) && isfinite(sigma) && sigma > 0)
  khat = 0;
  return;
end

% Expected order statistics of GPD excesses, then right-truncate
probs = ((1:M) - 0.5) / M;
if abs(khat) < 1e-10
  qq = -sigma * log1p(-probs);
else
  qq = (sigma / khat) * ((1 - probs) .^ (-khat) - 1);
end
qq = min(qq, S^(3/4));
newTail = cutoff + qq;
logW(ord(S - M + 1:S)) = newTail;
end
