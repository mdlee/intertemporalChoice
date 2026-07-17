function pmfPad = postPredictiveCountPmfFromTheta(Th, pairRow, nPairs, nTpRow, nPredSamples)
%POSTPREDICTIVECOUNTPMFFROMTHETA  Kruschke-style count PMF per problem (padded).
%
%   Th: nDraws x nTrials matrix of Bernoulli probabilities for one participant.

nkMax = max(nTpRow);
if nkMax <= 0 || isempty(Th)
  pmfPad = [];
  return;
end

nT = size(Th, 2);
if size(Th, 1) > nPredSamples && nPredSamples > 0
  idx = randperm(size(Th, 1), nPredSamples);
  Th = Th(idx, :);
end

Th = min(max(Th, 0), 1);
useS = size(Th, 1);
pmfPad = nan(nPairs, nkMax + 1);

cntSim = nan(useS, nPairs);
for s = 1:useS
  u = rand(1, nT);
  ySim = u < Th(s, :);
  for p = 1:nPairs
    m = pairRow == p;
    if any(m)
      cntSim(s, p) = sum(ySim(m));
    end
  end
end

for p = 1:nPairs
  nk = nTpRow(p);
  if nk <= 0
    continue;
  end
  col = cntSim(:, p);
  col = col(isfinite(col));
  if isempty(col)
    continue;
  end
  for k = 0:nk
    pmfPad(p, k + 1) = mean(col == k);
  end
end

end
