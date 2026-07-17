function Th = hierarchicalExecutionThetaDraws(chains, data, modelName, pp)
%HIERARCHICALEXECUTIONTHETADRAWS  nDraws-by-nTrials theta for one participant.
%
%   Computes decision probabilities from hierarchical CODA
%   (param_1, param_2, ...) without monitored theta_* fields.
%   Models ending in _entrop use Grünwald entropification (w per participant).

rLL = data.rLL(pp, :);
rSS = data.rSS(pp, :);
tLL = data.tLL(pp, :);
tSS = data.tSS(pp, :);

useEntrop = endsWith(modelName, '_entrop');
if useEntrop
  modelName = extractBefore(modelName, '_entrop');
end
% Prior-robustness stems: ...MuPrecHalf_entrop / ...MuPrecDouble_entrop
modelName = strrep(modelName, 'MuPrecHalf', '');
modelName = strrep(modelName, 'MuPrecDouble', '');

switch modelName
  case 'exponentialExecutionHierarchical'
    kappa = participantDraws(chains, 'kappa', pp);
    vLL = rLL .* exp(-kappa .* tLL);
    vSS = rSS .* exp(-kappa .* tSS);
    if useEntrop
      Th = entropThetaFromLatent(vLL - vSS, participantDraws(chains, 'w', pp));
    else
      Th = executionThetaFromDiff(vLL - vSS, participantDraws(chains, 'epsilon', pp));
    end

  case 'hyperbolicExecutionHierarchical'
    kappa = participantDraws(chains, 'kappa', pp);
    vLL = rLL ./ (1 + kappa .* tLL);
    vSS = rSS ./ (1 + kappa .* tSS);
    if useEntrop
      Th = entropThetaFromLatent(vLL - vSS, participantDraws(chains, 'w', pp));
    else
      Th = executionThetaFromDiff(vLL - vSS, participantDraws(chains, 'epsilon', pp));
    end

  case 'hyperboloidExecutionHierarchical'
    kappa = participantDraws(chains, 'kappa', pp);
    tau = participantDraws(chains, 'tau', pp);
    vLL = rLL ./ (1 + kappa .* tLL) .^ tau;
    vSS = rSS ./ (1 + kappa .* tSS) .^ tau;
    if useEntrop
      Th = entropThetaFromLatent(vLL - vSS, participantDraws(chains, 'w', pp));
    else
      Th = executionThetaFromDiff(vLL - vSS, participantDraws(chains, 'epsilon', pp));
    end

  case 'proportionalDifferencesExecutionHierarchical'
    delta = participantDraws(chains, 'delta', pp);
    dr = (rLL - rSS) ./ rLL;
    dt = (tLL - tSS) ./ tLL;
    if useEntrop
      Th = entropThetaFromLatent(dr - dt - delta, participantDraws(chains, 'w', pp));
    else
      Th = executionThetaFromDiff(dr - dt - delta, participantDraws(chains, 'epsilon', pp));
    end

  case 'directDifferencesExecutionHierarchical'
    delta = participantDraws(chains, 'delta', pp);
    omega = participantDraws(chains, 'omega', pp);
    dr = omega .* (rLL - rSS);
    dt = (1 - omega) .* (tLL - tSS);
    if useEntrop
      Th = entropThetaFromLatent(dr - dt - delta, participantDraws(chains, 'w', pp));
    else
      Th = executionThetaFromDiff(dr - dt - delta, participantDraws(chains, 'epsilon', pp));
    end

  case 'tradeoffExecutionHierarchical'
    gamma = participantDraws(chains, 'gamma', pp);
    tau = participantDraws(chains, 'tau', pp);
    kappa = participantDraws(chains, 'kappa', pp);
    vartheta = participantDraws(chains, 'vartheta', pp);
    if useEntrop
      Th = entropTradeoffThetaDraws(gamma, tau, kappa, vartheta, ...
        participantDraws(chains, 'w', pp), rLL, rSS, tLL, tSS, false, []);
    else
      Th = tradeoffThetaDraws(gamma, tau, kappa, vartheta, ...
        participantDraws(chains, 'epsilon', pp), rLL, rSS, tLL, tSS, false, []);
    end

  case 'unifiedTradeoffExecutionHierarchical'
    gamma = participantDraws(chains, 'gamma', pp);
    tau = participantDraws(chains, 'tau', pp);
    kappa = participantDraws(chains, 'kappa', pp);
    vartheta = participantDraws(chains, 'vartheta', pp);
    eta = participantDraws(chains, 'eta', pp);
    if useEntrop
      Th = entropTradeoffThetaDraws(gamma, tau, kappa, vartheta, ...
        participantDraws(chains, 'w', pp), rLL, rSS, tLL, tSS, true, eta);
    else
      Th = tradeoffThetaDraws(gamma, tau, kappa, vartheta, ...
        participantDraws(chains, 'epsilon', pp), rLL, rSS, tLL, tSS, true, eta);
    end

  case 'itchExecutionHierarchical'
    beta0 = participantDraws(chains, 'beta0', pp);
    betaRA = participantDraws(chains, 'betaRA', pp);
    betaRR = participantDraws(chains, 'betaRR', pp);
    betaTA = participantDraws(chains, 'betaTA', pp);
    betaTR = participantDraws(chains, 'betaTR', pp);
    drA = betaRA .* (rLL - rSS);
    drR = betaRR .* ((rLL - rSS) ./ (0.5 * (rLL + rSS)));
    dtA = betaTA .* (tLL - tSS);
    dtR = betaTR .* ((tLL - tSS) ./ (0.5 * (tLL + tSS)));
    score = beta0 + drA + drR + dtA + dtR;
    if useEntrop
      Th = entropThetaFromLatent(score, participantDraws(chains, 'w', pp));
    else
      Th = executionThetaFromDiff(score, participantDraws(chains, 'epsilon', pp));
    end

  case 'exponentialHierarchical'
    kappa = participantDraws(chains, 'kappa', pp);
    vLL = rLL .* exp(-kappa .* tLL);
    vSS = rSS .* exp(-kappa .* tSS);
    if useEntrop
      Th = entropThetaFromLatent(vLL - vSS, participantDraws(chains, 'w', pp));
    else
      Th = probitThetaFromLatent(vLL - vSS, participantDraws(chains, 'sigma', pp));
    end

  case 'hyperbolicHierarchical'
    kappa = participantDraws(chains, 'kappa', pp);
    vLL = rLL ./ (1 + kappa .* tLL);
    vSS = rSS ./ (1 + kappa .* tSS);
    if useEntrop
      Th = entropThetaFromLatent(vLL - vSS, participantDraws(chains, 'w', pp));
    else
      Th = probitThetaFromLatent(vLL - vSS, participantDraws(chains, 'sigma', pp));
    end

  case 'hyperboloidHierarchical'
    kappa = participantDraws(chains, 'kappa', pp);
    tau = participantDraws(chains, 'tau', pp);
    vLL = rLL ./ (1 + kappa .* tLL) .^ tau;
    vSS = rSS ./ (1 + kappa .* tSS) .^ tau;
    if useEntrop
      Th = entropThetaFromLatent(vLL - vSS, participantDraws(chains, 'w', pp));
    else
      Th = probitThetaFromLatent(vLL - vSS, participantDraws(chains, 'sigma', pp));
    end

  case 'proportionalDifferencesHierarchical'
    delta = participantDraws(chains, 'delta', pp);
    dr = (rLL - rSS) ./ rLL;
    dt = (tLL - tSS) ./ tLL;
    if useEntrop
      Th = entropThetaFromLatent(dr - dt - delta, participantDraws(chains, 'w', pp));
    else
      Th = probitThetaFromLatent(dr - dt - delta, participantDraws(chains, 'sigma', pp));
    end

  case 'directDifferencesHierarchical'
    delta = participantDraws(chains, 'delta', pp);
    omega = participantDraws(chains, 'omega', pp);
    dr = omega .* (rLL - rSS);
    dt = (1 - omega) .* (tLL - tSS);
    if useEntrop
      Th = entropThetaFromLatent(dr - dt - delta, participantDraws(chains, 'w', pp));
    else
      Th = probitThetaFromLatent(dr - dt - delta, participantDraws(chains, 'sigma', pp));
    end

  case 'tradeoffHierarchical'
    gamma = participantDraws(chains, 'gamma', pp);
    tau = participantDraws(chains, 'tau', pp);
    kappa = participantDraws(chains, 'kappa', pp);
    vartheta = participantDraws(chains, 'vartheta', pp);
    if useEntrop
      Th = entropTradeoffThetaDraws(gamma, tau, kappa, vartheta, ...
        participantDraws(chains, 'w', pp), rLL, rSS, tLL, tSS, false, []);
    else
      Th = probitTradeoffThetaDraws(gamma, tau, kappa, vartheta, ...
        participantDraws(chains, 'sigma', pp), rLL, rSS, tLL, tSS, false, []);
    end

  case 'unifiedTradeoffHierarchical'
    gamma = participantDraws(chains, 'gamma', pp);
    tau = participantDraws(chains, 'tau', pp);
    kappa = participantDraws(chains, 'kappa', pp);
    vartheta = participantDraws(chains, 'vartheta', pp);
    eta = participantDraws(chains, 'eta', pp);
    if useEntrop
      Th = entropTradeoffThetaDraws(gamma, tau, kappa, vartheta, ...
        participantDraws(chains, 'w', pp), rLL, rSS, tLL, tSS, true, eta);
    else
      Th = probitTradeoffThetaDraws(gamma, tau, kappa, vartheta, ...
        participantDraws(chains, 'sigma', pp), rLL, rSS, tLL, tSS, true, eta);
    end

  case 'itchHierarchical'
    beta0 = participantDraws(chains, 'beta0', pp);
    betaRA = participantDraws(chains, 'betaRA', pp);
    betaRR = participantDraws(chains, 'betaRR', pp);
    betaTA = participantDraws(chains, 'betaTA', pp);
    betaTR = participantDraws(chains, 'betaTR', pp);
    drA = betaRA .* (rLL - rSS);
    drR = betaRR .* ((rLL - rSS) ./ (0.5 * (rLL + rSS)));
    dtA = betaTA .* (tLL - tSS);
    dtR = betaTR .* ((tLL - tSS) ./ (0.5 * (tLL + tSS)));
    score = beta0 + drA + drR + dtA + dtR;
    if useEntrop
      Th = entropThetaFromLatent(score, participantDraws(chains, 'w', pp));
    else
      Th = probitThetaFromLatent(score, participantDraws(chains, 'sigma', pp));
    end

  otherwise
    error('hierarchicalExecutionThetaDraws:unknownModel', ...
      'No theta mapping for model %s', modelName);
end

Th = min(max(Th, 0), 1);

end

function x = participantDraws(chains, paramName, pp)
mats = codaIndexedMatrices(chains, paramName);
x = mats{pp}(:);
end

function Th = executionThetaFromDiff(diffV, epsilon)
% diffV: nDraws x nTrials, epsilon: nDraws x 1
Th = epsilon + (1 - 2 * epsilon) .* double(diffV > 0);
end

function Th = probitThetaFromLatent(latent, sigma)
% latent: nDraws x nTrials, sigma: nDraws x 1
Th = normcdf(latent ./ max(sigma, 1e-12));
end

function Th = entropThetaFromLatent(latent, w)
% latent: nDraws x nTrials, w: nDraws x 1
nD = size(latent, 1);
nT = size(latent, 2);
Th = nan(nD, nT);
pLL = 1 ./ (1 + exp(-w));
pSS = 1 ./ (1 + exp(w));
for s = 1:nD
  pick = latent(s, :) > 0;
  Th(s, pick) = pLL(s);
  Th(s, ~pick) = pSS(s);
end
end

function Th = probitTradeoffThetaDraws(gamma, tau, kappa, vartheta, sigma, rLL, rSS, tLL, tSS, unified, eta)
nD = numel(gamma);
nT = numel(rLL);
Th = nan(nD, nT);
for s = 1:nD
  g = gamma(s); tv = tau(s); k = kappa(s); vt = vartheta(s); sg = sigma(s);
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
  Th(s, :) = normcdf((qV - qT) / max(sg, 1e-12));
end
end

function Th = entropTradeoffThetaDraws(gamma, tau, kappa, vartheta, w, rLL, rSS, tLL, tSS, unified, eta)
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
  latentRow = qV - qT;
  pick = latentRow > 0;
  pLL = 1 / (1 + exp(-w(s)));
  pSS = 1 / (1 + exp(w(s)));
  Th(s, pick) = pLL;
  Th(s, ~pick) = pSS;
end
end

function Th = tradeoffThetaDraws(gamma, tau, kappa, vartheta, epsilon, rLL, rSS, tLL, tSS, unified, eta)
nD = numel(gamma);
nT = numel(rLL);
Th = nan(nD, nT);
for s = 1:nD
  g = gamma(s); tv = tau(s); k = kappa(s); vt = vartheta(s); e = epsilon(s);
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
  diffRow = qV - qT;
  Th(s, :) = e + (1 - 2 * e) * double(diffRow > 0);
end
end
