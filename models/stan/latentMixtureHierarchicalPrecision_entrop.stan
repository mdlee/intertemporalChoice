// 11-component latent mixture (manuscript): 8 hierarchical cognitive models
// + Guess / LL / SS contaminants. One z per participant (marginalized).
// Cognitive P(LL) uses Grünwald entropification via w.
// Robustness: scale group-mean prior precision with muPrecScale (1, 0.5, 2).

functions {
  real entrop_logit(real advantage, real w) {
    return w * (2 * step(advantage) - 1);
  }

  vector mix_trial_loglik(int y, real rLL, real rSS, real tLL, real tSS,
      real kappaEX, real kappaHC, real kappaHY, real tauHY,
      real deltaPD, real deltaDD, real omegaDD,
      real gammaTR, real tauTR, real kappaTR, real varthetaTR,
      real gammaUT, real tauUT, real kappaUT, real varthetaUT, real etaUT,
      real beta0IT, real betaRAIT, real betaRRIT, real betaTAIT, real betaTRIT,
      real w, real alphaLL, real alphaSS) {
    vector[11] lp;
    real adv;
    real vLL;
    real vSS;
    real wLL;
    real wSS;
    real qV;
    real dW;
    real qT;
    real rNet;

    vLL = rLL * exp(-kappaEX * tLL);
    vSS = rSS * exp(-kappaEX * tSS);
    lp[1] = bernoulli_logit_lpmf(y | entrop_logit(vLL - vSS, w));

    vLL = rLL / (1 + kappaHC * tLL);
    vSS = rSS / (1 + kappaHC * tSS);
    lp[2] = bernoulli_logit_lpmf(y | entrop_logit(vLL - vSS, w));

    vLL = rLL / pow(1 + kappaHY * tLL, tauHY);
    vSS = rSS / pow(1 + kappaHY * tSS, tauHY);
    lp[3] = bernoulli_logit_lpmf(y | entrop_logit(vLL - vSS, w));

    adv = ((rLL - rSS) / rLL) - ((tLL - tSS) / tLL) - deltaPD;
    lp[4] = bernoulli_logit_lpmf(y | entrop_logit(adv, w));

    adv = omegaDD * (rLL - rSS) - (1 - omegaDD) * (tLL - tSS) - deltaDD;
    lp[5] = bernoulli_logit_lpmf(y | entrop_logit(adv, w));

    vLL = (1 / gammaTR) * log(1 + gammaTR * rLL);
    vSS = (1 / gammaTR) * log(1 + gammaTR * rSS);
    wLL = (1 / tauTR) * log(1 + tauTR * tLL);
    wSS = (1 / tauTR) * log(1 + tauTR * tSS);
    qV = vLL - vSS;
    dW = fmax(wLL - wSS, 0);
    qT = kappaTR * log(1 + pow(dW / varthetaTR, varthetaTR));
    lp[6] = bernoulli_logit_lpmf(y | entrop_logit(qV - qT, w));

    rNet = fmax(rLL - etaUT, 1e-6);
    vLL = (1 / gammaUT) * log(1 + gammaUT * rNet);
    vSS = (1 / gammaUT) * log(1 + gammaUT * rSS);
    wLL = (1 / tauUT) * log(1 + tauUT * tLL);
    wSS = (1 / tauUT) * log(1 + tauUT * tSS);
    qV = vLL - vSS;
    dW = fmax(wLL - wSS, 0);
    qT = kappaUT * log(1 + pow(dW / varthetaUT, varthetaUT));
    lp[7] = bernoulli_logit_lpmf(y | entrop_logit(qV - qT, w));

    adv = beta0IT
        + betaRAIT * (rLL - rSS)
        + betaRRIT * ((rLL - rSS) / (0.5 * (rLL + rSS)))
        + betaTAIT * (tLL - tSS)
        + betaTRIT * ((tLL - tSS) / (0.5 * (tLL + tSS)));
    lp[8] = bernoulli_logit_lpmf(y | entrop_logit(adv, w));

    lp[9] = bernoulli_lpmf(y | 0.5);
    lp[10] = bernoulli_lpmf(y | alphaLL);
    lp[11] = bernoulli_lpmf(y | 1 - alphaSS);
    return lp;
  }
}

data {
  int<lower=1> nParticipants;
  int<lower=1> nTrials;
  int<lower=1> nModels;
  array[nParticipants, nTrials] int<lower=0, upper=1> decision;
  matrix[nParticipants, nTrials] rLL;
  matrix[nParticipants, nTrials] rSS;
  matrix[nParticipants, nTrials] tLL;
  matrix[nParticipants, nTrials] tSS;
  real<lower=0> muPrecScale;
}

parameters {
  real<lower=1e-4, upper=100> muKappaEX;
  real<lower=0.05, upper=25> tauKappaEX;
  real<lower=1e-4, upper=100> muKappaHC;
  real<lower=0.05, upper=25> tauKappaHC;
  real<lower=1e-4, upper=100> muKappaHY;
  real<lower=0.05, upper=25> tauKappaHY;
  real<lower=1e-4, upper=100> muTauHY;
  real<lower=0.05, upper=25> tauTauHY;
  real<lower=-2, upper=2> muDeltaPD;
  real<lower=0.05, upper=25> tauDeltaPD;
  real<lower=-2000, upper=2000> muDeltaDD;
  real<lower=0.05, upper=25> tauDeltaDD;
  real<lower=1e-6, upper=1 - 1e-6> muOmegaDD;
  real<lower=0.05, upper=25> tauOmegaDD;
  real<lower=1e-4, upper=100> muGammaTR;
  real<lower=0.05, upper=25> tauGammaTR;
  real<lower=1e-4, upper=100> muTauTR;
  real<lower=0.05, upper=25> tauTauTR;
  real<lower=1e-4, upper=100> muKappaTR;
  real<lower=0.05, upper=25> tauKappaTR;
  real<lower=1e-4, upper=100> muVarthetaTR;
  real<lower=0.05, upper=25> tauVarthetaTR;
  real<lower=1e-4, upper=100> muGammaUT;
  real<lower=0.05, upper=25> tauGammaUT;
  real<lower=1e-4, upper=100> muTauUT;
  real<lower=0.05, upper=25> tauTauUT;
  real<lower=1e-4, upper=100> muKappaUT;
  real<lower=0.05, upper=25> tauKappaUT;
  real<lower=1e-4, upper=100> muVarthetaUT;
  real<lower=0.05, upper=25> tauVarthetaUT;
  real<lower=0, upper=500> muEtaUT;
  real<lower=0.05, upper=25> tauEtaUT;
  real<lower=-2, upper=2> muBeta0IT;
  real<lower=0.05, upper=25> tauBeta0IT;
  real<lower=0, upper=2> muBetaRAIT;
  real<lower=0.05, upper=25> tauBetaRAIT;
  real<lower=0, upper=2> muBetaRRIT;
  real<lower=0.05, upper=25> tauBetaRRIT;
  real<lower=-2, upper=0> muBetaTAIT;
  real<lower=0.05, upper=25> tauBetaTAIT;
  real<lower=-2, upper=0> muBetaTRIT;
  real<lower=0.05, upper=25> tauBetaTRIT;
  real<lower=0, upper=20> muW;
  real<lower=0.05, upper=25> tauW;

  vector<lower=1e-4, upper=100>[nParticipants] kappaEX;
  vector<lower=1e-4, upper=100>[nParticipants] kappaHC;
  vector<lower=1e-4, upper=100>[nParticipants] kappaHY;
  vector<lower=1e-4, upper=100>[nParticipants] tauHY;
  vector<lower=-2, upper=2>[nParticipants] deltaPD;
  vector<lower=-2000, upper=2000>[nParticipants] deltaDD;
  vector<lower=1e-6, upper=1 - 1e-6>[nParticipants] omegaDD;
  vector<lower=1e-4, upper=100>[nParticipants] gammaTR;
  vector<lower=1e-4, upper=100>[nParticipants] tauTR;
  vector<lower=1e-4, upper=100>[nParticipants] kappaTR;
  vector<lower=1e-4, upper=100>[nParticipants] varthetaTR;
  vector<lower=1e-4, upper=100>[nParticipants] gammaUT;
  vector<lower=1e-4, upper=100>[nParticipants] tauUT;
  vector<lower=1e-4, upper=100>[nParticipants] kappaUT;
  vector<lower=1e-4, upper=100>[nParticipants] varthetaUT;
  vector<lower=0, upper=500>[nParticipants] etaUT;
  vector<lower=-2, upper=2>[nParticipants] beta0IT;
  vector<lower=0, upper=2>[nParticipants] betaRAIT;
  vector<lower=0, upper=2>[nParticipants] betaRRIT;
  vector<lower=-2, upper=0>[nParticipants] betaTAIT;
  vector<lower=-2, upper=0>[nParticipants] betaTRIT;
  vector<lower=0, upper=20>[nParticipants] w;
  vector<lower=1e-6, upper=1 - 1e-6>[nParticipants] alphaLL;
  vector<lower=1e-6, upper=1 - 1e-6>[nParticipants] alphaSS;
}

model {
  muKappaEX ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauKappaEX ~ gamma(1, 1);
  muKappaHC ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauKappaHC ~ gamma(1, 1);
  muKappaHY ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauKappaHY ~ gamma(1, 1);
  muTauHY ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauTauHY ~ gamma(1, 1);
  muDeltaPD ~ normal(0, inv_sqrt(4 * muPrecScale));
  tauDeltaPD ~ gamma(1, 1);
  muDeltaDD ~ normal(0, inv_sqrt(1e-6 * muPrecScale));
  tauDeltaDD ~ gamma(1, 1);
  muOmegaDD ~ normal(0.5, inv_sqrt(4 * muPrecScale));
  tauOmegaDD ~ gamma(1, 1);
  muGammaTR ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauGammaTR ~ gamma(1, 1);
  muTauTR ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauTauTR ~ gamma(1, 1);
  muKappaTR ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauKappaTR ~ gamma(1, 1);
  muVarthetaTR ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauVarthetaTR ~ gamma(1, 1);
  muGammaUT ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauGammaUT ~ gamma(1, 1);
  muTauUT ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauTauUT ~ gamma(1, 1);
  muKappaUT ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauKappaUT ~ gamma(1, 1);
  muVarthetaUT ~ normal(1, inv_sqrt(4 * muPrecScale));
  tauVarthetaUT ~ gamma(1, 1);
  muEtaUT ~ normal(0, inv_sqrt(4 * muPrecScale));
  tauEtaUT ~ gamma(1, 1);
  muBeta0IT ~ normal(0, inv_sqrt(4 * muPrecScale));
  tauBeta0IT ~ gamma(1, 1);
  muBetaRAIT ~ normal(0, inv_sqrt(4 * muPrecScale));
  tauBetaRAIT ~ gamma(1, 1);
  muBetaRRIT ~ normal(0, inv_sqrt(4 * muPrecScale));
  tauBetaRRIT ~ gamma(1, 1);
  muBetaTAIT ~ normal(0, inv_sqrt(4 * muPrecScale));
  tauBetaTAIT ~ gamma(1, 1);
  muBetaTRIT ~ normal(0, inv_sqrt(4 * muPrecScale));
  tauBetaTRIT ~ gamma(1, 1);
  muW ~ normal(2, inv_sqrt(1 * muPrecScale));
  tauW ~ gamma(1, 1);

  kappaEX ~ normal(muKappaEX, inv_sqrt(tauKappaEX));
  kappaHC ~ normal(muKappaHC, inv_sqrt(tauKappaHC));
  kappaHY ~ normal(muKappaHY, inv_sqrt(tauKappaHY));
  tauHY ~ normal(muTauHY, inv_sqrt(tauTauHY));
  deltaPD ~ normal(muDeltaPD, inv_sqrt(tauDeltaPD));
  deltaDD ~ normal(muDeltaDD, inv_sqrt(tauDeltaDD));
  omegaDD ~ normal(muOmegaDD, inv_sqrt(tauOmegaDD));
  gammaTR ~ normal(muGammaTR, inv_sqrt(tauGammaTR));
  tauTR ~ normal(muTauTR, inv_sqrt(tauTauTR));
  kappaTR ~ normal(muKappaTR, inv_sqrt(tauKappaTR));
  varthetaTR ~ normal(muVarthetaTR, inv_sqrt(tauVarthetaTR));
  gammaUT ~ normal(muGammaUT, inv_sqrt(tauGammaUT));
  tauUT ~ normal(muTauUT, inv_sqrt(tauTauUT));
  kappaUT ~ normal(muKappaUT, inv_sqrt(tauKappaUT));
  varthetaUT ~ normal(muVarthetaUT, inv_sqrt(tauVarthetaUT));
  etaUT ~ normal(muEtaUT, inv_sqrt(tauEtaUT));
  beta0IT ~ normal(muBeta0IT, inv_sqrt(tauBeta0IT));
  betaRAIT ~ normal(muBetaRAIT, inv_sqrt(tauBetaRAIT));
  betaRRIT ~ normal(muBetaRRIT, inv_sqrt(tauBetaRRIT));
  betaTAIT ~ normal(muBetaTAIT, inv_sqrt(tauBetaTAIT));
  betaTRIT ~ normal(muBetaTRIT, inv_sqrt(tauBetaTRIT));
  w ~ normal(muW, inv_sqrt(tauW));
  alphaLL ~ beta(10, 1);
  alphaSS ~ beta(10, 1);

  for (i in 1:nParticipants) {
    vector[11] lp = rep_vector(-log(nModels), 11);
    for (j in 1:nTrials) {
      lp += mix_trial_loglik(decision[i, j], rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j],
          kappaEX[i], kappaHC[i], kappaHY[i], tauHY[i],
          deltaPD[i], deltaDD[i], omegaDD[i],
          gammaTR[i], tauTR[i], kappaTR[i], varthetaTR[i],
          gammaUT[i], tauUT[i], kappaUT[i], varthetaUT[i], etaUT[i],
          beta0IT[i], betaRAIT[i], betaRRIT[i], betaTAIT[i], betaTRIT[i],
          w[i], alphaLL[i], alphaSS[i]);
    }
    target += log_sum_exp(lp);
  }
}

generated quantities {
  array[nParticipants] int<lower=1, upper=11> z;
  for (i in 1:nParticipants) {
    vector[11] lp = rep_vector(-log(nModels), 11);
    for (j in 1:nTrials) {
      lp += mix_trial_loglik(decision[i, j], rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j],
          kappaEX[i], kappaHC[i], kappaHY[i], tauHY[i],
          deltaPD[i], deltaDD[i], omegaDD[i],
          gammaTR[i], tauTR[i], kappaTR[i], varthetaTR[i],
          gammaUT[i], tauUT[i], kappaUT[i], varthetaUT[i], etaUT[i],
          beta0IT[i], betaRAIT[i], betaRRIT[i], betaTAIT[i], betaTRIT[i],
          w[i], alphaLL[i], alphaSS[i]);
    }
    z[i] = categorical_logit_rng(lp);
  }
}
