function ll = latentMixturePointwiseLogLik(chains, data)
%LATENTMIXTUREPOINTWISELOGLIK  nObs-by-nDraws log Bern likelihood (one chain).
%
%   CHAINS fields are nDraws-by-1 (single chain). Observations are stacked as
%   participant-major: obs = (pp-1)*nTrials + trial.
%   Matches latentMixtureHierarchicalPrecision_entrop_jags.txt.

nP = data.nParticipants;
nT = data.nTrials;
nS = size(chains.z_1, 1);
ll = nan(nP * nT, nS);

for pp = 1:nP
  z = max(1, min(11, round(chains.(sprintf('z_%d', pp))(:))));
  w = chains.(sprintf('w_%d', pp)); w = w(:);
  pLL = 1 ./ (1 + exp(-w));
  pSS = 1 ./ (1 + exp(w));

  rLL = data.rLL(pp, :);
  rSS = data.rSS(pp, :);
  tLL = data.tLL(pp, :);
  tSS = data.tSS(pp, :);
  y = data.decision(pp, :);

  kappaEX = chains.(sprintf('kappaEX_%d', pp)); kappaEX = kappaEX(:);
  kappaHC = chains.(sprintf('kappaHC_%d', pp)); kappaHC = kappaHC(:);
  kappaHY = chains.(sprintf('kappaHY_%d', pp)); kappaHY = kappaHY(:);
  tauHY = chains.(sprintf('tauHY_%d', pp)); tauHY = tauHY(:);
  deltaPD = chains.(sprintf('deltaPD_%d', pp)); deltaPD = deltaPD(:);
  deltaDD = chains.(sprintf('deltaDD_%d', pp)); deltaDD = deltaDD(:);
  omegaDD = chains.(sprintf('omegaDD_%d', pp)); omegaDD = omegaDD(:);
  gammaTR = chains.(sprintf('gammaTR_%d', pp)); gammaTR = gammaTR(:);
  tauTR = chains.(sprintf('tauTR_%d', pp)); tauTR = tauTR(:);
  kappaTR = chains.(sprintf('kappaTR_%d', pp)); kappaTR = kappaTR(:);
  varthetaTR = chains.(sprintf('varthetaTR_%d', pp)); varthetaTR = varthetaTR(:);
  gammaUT = chains.(sprintf('gammaUT_%d', pp)); gammaUT = gammaUT(:);
  tauUT = chains.(sprintf('tauUT_%d', pp)); tauUT = tauUT(:);
  kappaUT = chains.(sprintf('kappaUT_%d', pp)); kappaUT = kappaUT(:);
  varthetaUT = chains.(sprintf('varthetaUT_%d', pp)); varthetaUT = varthetaUT(:);
  etaUT = chains.(sprintf('etaUT_%d', pp)); etaUT = etaUT(:);
  beta0IT = chains.(sprintf('beta0IT_%d', pp)); beta0IT = beta0IT(:);
  betaRAIT = chains.(sprintf('betaRAIT_%d', pp)); betaRAIT = betaRAIT(:);
  betaRRIT = chains.(sprintf('betaRRIT_%d', pp)); betaRRIT = betaRRIT(:);
  betaTAIT = chains.(sprintf('betaTAIT_%d', pp)); betaTAIT = betaTAIT(:);
  betaTRIT = chains.(sprintf('betaTRIT_%d', pp)); betaTRIT = betaTRIT(:);
  alphaLL = chains.(sprintf('alphaLL_%d', pp)); alphaLL = alphaLL(:);
  alphaSS = chains.(sprintf('alphaSS_%d', pp)); alphaSS = alphaSS(:);

  % nS × nT latent advantages for each cognitive model
  L = zeros(nS, nT, 8);
  L(:, :, 1) = rLL .* exp(-kappaEX .* tLL) - rSS .* exp(-kappaEX .* tSS);
  L(:, :, 2) = rLL ./ (1 + kappaHC .* tLL) - rSS ./ (1 + kappaHC .* tSS);
  L(:, :, 3) = rLL ./ (1 + kappaHY .* tLL) .^ tauHY ...
             - rSS ./ (1 + kappaHY .* tSS) .^ tauHY;
  L(:, :, 4) = (rLL - rSS) ./ rLL - (tLL - tSS) ./ tLL - deltaPD;
  L(:, :, 5) = omegaDD .* (rLL - rSS) - (1 - omegaDD) .* (tLL - tSS) - deltaDD;

  vLLtr = (1 ./ gammaTR) .* log(1 + gammaTR .* rLL);
  vSStr = (1 ./ gammaTR) .* log(1 + gammaTR .* rSS);
  wLLtr = (1 ./ tauTR) .* log(1 + tauTR .* tLL);
  wSStr = (1 ./ tauTR) .* log(1 + tauTR .* tSS);
  L(:, :, 6) = (vLLtr - vSStr) - kappaTR .* log(1 + max(((wLLtr - wSStr) ./ varthetaTR) .^ varthetaTR, 0));

  vLLut = (1 ./ gammaUT) .* log(1 + gammaUT .* max(rLL - etaUT, 1e-9));
  vSSut = (1 ./ gammaUT) .* log(1 + gammaUT .* rSS);
  wLLut = (1 ./ tauUT) .* log(1 + tauUT .* tLL);
  wSSut = (1 ./ tauUT) .* log(1 + tauUT .* tSS);
  L(:, :, 7) = (vLLut - vSSut) - kappaUT .* log(1 + max(((wLLut - wSSut) ./ varthetaUT) .^ varthetaUT, 0));

  L(:, :, 8) = beta0IT ...
    + betaRAIT .* (rLL - rSS) ...
    + betaRRIT .* ((rLL - rSS) ./ (0.5 * (rLL + rSS))) ...
    + betaTAIT .* (tLL - tSS) ...
    + betaTRIT .* ((tLL - tSS) ./ (0.5 * (tLL + tSS)));

  Th = zeros(nS, nT);
  for s = 1:nS
    m = z(s);
    if m <= 8
      pick = L(s, :, m) > 0;
      row = zeros(1, nT);
      row(pick) = pLL(s);
      row(~pick) = pSS(s);
      Th(s, :) = row;
    elseif m == 9
      Th(s, :) = 0.5;
    elseif m == 10
      Th(s, :) = alphaLL(s);
    else
      Th(s, :) = 1 - alphaSS(s);
    end
  end
  Th = min(max(Th, 1e-12), 1 - 1e-12);

  rows = (pp - 1) * nT + (1:nT);
  for j = 1:nT
    if y(j) == 1
      ll(rows(j), :) = log(Th(:, j))';
    else
      ll(rows(j), :) = log(1 - Th(:, j))';
    end
  end
end
end
