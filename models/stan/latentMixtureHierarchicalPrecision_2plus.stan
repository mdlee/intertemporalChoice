#include "common_functions.stan"

data {
  int<lower=1> nParticipants;
  int<lower=1> nTrials;
  int<lower=1> nModels;
  matrix[nParticipants, nTrials] rLL;
  matrix[nParticipants, nTrials] rSS;
  matrix[nParticipants, nTrials] tLL;
  matrix[nParticipants, nTrials] tSS;
  array[nParticipants, nTrials] int<lower=0, upper=1> decision;
}
parameters {
  real<lower=1e-4, upper=100> mu_kappaEX;
  real<lower=0.05, upper=25> tau_kappaEX;
  real<lower=0, upper=0.5> mu_epsilonEX;
  real<lower=0.05, upper=25> tau_epsilonEX;
  real<lower=1e-4, upper=100> mu_kappaHC;
  real<lower=0.05, upper=25> tau_kappaHC;
  real<lower=0, upper=0.5> mu_epsilonHC;
  real<lower=0.05, upper=25> tau_epsilonHC;
  real<lower=1e-4, upper=100> mu_kappaHY;
  real<lower=0.05, upper=25> tau_kappaHY;
  real<lower=1e-4, upper=100> mu_tauHY;
  real<lower=0.05, upper=25> tau_tauHY;
  real<lower=0, upper=0.5> mu_epsilonHY;
  real<lower=0.05, upper=25> tau_epsilonHY;
  real<lower=-2, upper=2> mu_deltaPD;
  real<lower=0.05, upper=25> tau_deltaPD;
  real<lower=0, upper=0.5> mu_epsilonPD;
  real<lower=0.05, upper=25> tau_epsilonPD;
  real<lower=-2000, upper=2000> mu_deltaDD;
  real<lower=0.05, upper=25> tau_deltaDD;
  real<lower=0.01, upper=0.99> mu_omegaDD;
  real<lower=0.05, upper=25> tau_omegaDD;
  real<lower=0, upper=0.5> mu_epsilonDD;
  real<lower=0.05, upper=25> tau_epsilonDD;
  real<lower=1e-4, upper=100> mu_gammaTR;
  real<lower=0.05, upper=25> tau_gammaTR;
  real<lower=1e-4, upper=100> mu_tauTR;
  real<lower=0.05, upper=25> tau_tauTR;
  real<lower=1e-4, upper=100> mu_kappaTR;
  real<lower=0.05, upper=25> tau_kappaTR;
  real<lower=1e-4, upper=100> mu_varthetaTR;
  real<lower=0.05, upper=25> tau_varthetaTR;
  real<lower=0, upper=0.5> mu_epsilonTR;
  real<lower=0.05, upper=25> tau_epsilonTR;
  real<lower=1e-4, upper=100> mu_gammaUT;
  real<lower=0.05, upper=25> tau_gammaUT;
  real<lower=1e-4, upper=100> mu_tauUT;
  real<lower=0.05, upper=25> tau_tauUT;
  real<lower=1e-4, upper=100> mu_kappaUT;
  real<lower=0.05, upper=25> tau_kappaUT;
  real<lower=1e-4, upper=100> mu_varthetaUT;
  real<lower=0.05, upper=25> tau_varthetaUT;
  real<lower=0, upper=500> mu_etaUT;
  real<lower=0.05, upper=25> tau_etaUT;
  real<lower=0, upper=0.5> mu_epsilonUT;
  real<lower=0.05, upper=25> tau_epsilonUT;
  real<lower=-2, upper=2> mu_beta0IT;
  real<lower=0.05, upper=25> tau_beta0IT;
  real<lower=-2, upper=2> mu_betaRAIT;
  real<lower=0.05, upper=25> tau_betaRAIT;
  real<lower=-2, upper=2> mu_betaRRIT;
  real<lower=0.05, upper=25> tau_betaRRIT;
  real<lower=-2, upper=2> mu_betaTAIT;
  real<lower=0.05, upper=25> tau_betaTAIT;
  real<lower=-2, upper=2> mu_betaTRIT;
  real<lower=0.05, upper=25> tau_betaTRIT;
  real<lower=0, upper=0.5> mu_epsilonIT;
  real<lower=0.05, upper=25> tau_epsilonIT;
  real<lower=0.9, upper=1> mu_alphaLL;
  real<lower=0.05, upper=25> tau_alphaLL;
  real<lower=0.9, upper=1> mu_alphaSS;
  real<lower=0.05, upper=25> tau_alphaSS;
  vector<lower=1e-4, upper=100>[nParticipants] kappaEX;
  vector<lower=1e-4, upper=0.499>[nParticipants] epsilonEX;
  vector<lower=1e-4, upper=100>[nParticipants] kappaHC;
  vector<lower=1e-4, upper=0.499>[nParticipants] epsilonHC;
  vector<lower=1e-4, upper=100>[nParticipants] kappaHY;
  vector<lower=1e-4, upper=100>[nParticipants] tauHY;
  vector<lower=1e-4, upper=0.499>[nParticipants] epsilonHY;
  vector<lower=-2, upper=2>[nParticipants] deltaPD;
  vector<lower=1e-4, upper=0.499>[nParticipants] epsilonPD;
  vector<lower=-2000, upper=2000>[nParticipants] deltaDD;
  vector<lower=0.01, upper=0.99>[nParticipants] omegaDD;
  vector<lower=1e-4, upper=0.499>[nParticipants] epsilonDD;
  vector<lower=1e-4, upper=100>[nParticipants] gammaTR;
  vector<lower=1e-4, upper=100>[nParticipants] tauTR;
  vector<lower=1e-4, upper=100>[nParticipants] kappaTR;
  vector<lower=1e-4, upper=100>[nParticipants] varthetaTR;
  vector<lower=1e-4, upper=0.499>[nParticipants] epsilonTR;
  vector<lower=1e-4, upper=100>[nParticipants] gammaUT;
  vector<lower=1e-4, upper=100>[nParticipants] tauUT;
  vector<lower=1e-4, upper=100>[nParticipants] kappaUT;
  vector<lower=1e-4, upper=100>[nParticipants] varthetaUT;
  vector<lower=0, upper=500>[nParticipants] etaUT;
  vector<lower=1e-4, upper=0.499>[nParticipants] epsilonUT;
  vector<lower=-2, upper=2>[nParticipants] beta0IT;
  vector<lower=-2, upper=2>[nParticipants] betaRAIT;
  vector<lower=-2, upper=2>[nParticipants] betaRRIT;
  vector<lower=-2, upper=2>[nParticipants] betaTAIT;
  vector<lower=-2, upper=2>[nParticipants] betaTRIT;
  vector<lower=1e-4, upper=0.499>[nParticipants] epsilonIT;
  vector<lower=0.901, upper=0.999>[nParticipants] alphaLL;
  vector<lower=0.901, upper=0.999>[nParticipants] alphaSS;
}
transformed parameters {
  array[nParticipants, nTrials, nModels] real<lower=0, upper=1> theta;
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      theta[i, j, 1] = mix_ex_execution(
        rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j], kappaEX[i], epsilonEX[i]);
      theta[i, j, 2] = mix_hc_execution(
        rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j], kappaHC[i], epsilonHC[i]);
      theta[i, j, 3] = mix_hy_execution(
        rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j], kappaHY[i], tauHY[i], epsilonHY[i]);
      theta[i, j, 4] = mix_pd_execution(
        rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j], deltaPD[i], epsilonPD[i]);
      theta[i, j, 5] = mix_dd_execution(
        rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j], deltaDD[i], omegaDD[i], epsilonDD[i]);
      theta[i, j, 6] = mix_tr_execution(
        rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j],
        gammaTR[i], tauTR[i], kappaTR[i], varthetaTR[i], epsilonTR[i]);
      theta[i, j, 7] = mix_ut_execution(
        rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j],
        gammaUT[i], tauUT[i], kappaUT[i], varthetaUT[i], etaUT[i], epsilonUT[i]);
      theta[i, j, 8] = mix_it_execution(
        rLL[i, j], rSS[i, j], tLL[i, j], tSS[i, j],
        beta0IT[i], betaRAIT[i], betaRRIT[i], betaTAIT[i], betaTRIT[i], epsilonIT[i]);
      theta[i, j, 9] = 0.5;
      theta[i, j, 10] = alphaLL[i];
      theta[i, j, 11] = 1 - alphaSS[i];
    }
  }
}
model {
  real log_uniform;
  log_uniform = -log(nModels);
  mu_kappaEX ~ normal(1, 0.5);
  tau_kappaEX ~ gamma(1, 1);
  mu_epsilonEX ~ normal(0.05, 0.1);
  tau_epsilonEX ~ gamma(1, 1);
  kappaEX ~ normal(mu_kappaEX, 1 / sqrt(tau_kappaEX));
  epsilonEX ~ normal(mu_epsilonEX, 1 / sqrt(tau_epsilonEX));
  mu_kappaHC ~ normal(1, 0.5);
  tau_kappaHC ~ gamma(1, 1);
  mu_epsilonHC ~ normal(0.05, 0.1);
  tau_epsilonHC ~ gamma(1, 1);
  kappaHC ~ normal(mu_kappaHC, 1 / sqrt(tau_kappaHC));
  epsilonHC ~ normal(mu_epsilonHC, 1 / sqrt(tau_epsilonHC));
  mu_kappaHY ~ normal(1, 0.5);
  tau_kappaHY ~ gamma(1, 1);
  mu_tauHY ~ normal(1, 0.5);
  tau_tauHY ~ gamma(1, 1);
  mu_epsilonHY ~ normal(0.05, 0.1);
  tau_epsilonHY ~ gamma(1, 1);
  kappaHY ~ normal(mu_kappaHY, 1 / sqrt(tau_kappaHY));
  tauHY ~ normal(mu_tauHY, 1 / sqrt(tau_tauHY));
  epsilonHY ~ normal(mu_epsilonHY, 1 / sqrt(tau_epsilonHY));
  mu_deltaPD ~ normal(0, 0.5);
  tau_deltaPD ~ gamma(1, 1);
  mu_epsilonPD ~ normal(0.05, 0.1);
  tau_epsilonPD ~ gamma(1, 1);
  deltaPD ~ normal(mu_deltaPD, 1 / sqrt(tau_deltaPD));
  epsilonPD ~ normal(mu_epsilonPD, 1 / sqrt(tau_epsilonPD));
  mu_deltaDD ~ normal(0, 1000);
  tau_deltaDD ~ gamma(1, 1);
  mu_omegaDD ~ normal(0.5, 0.5);
  tau_omegaDD ~ gamma(1, 1);
  mu_epsilonDD ~ normal(0.05, 0.1);
  tau_epsilonDD ~ gamma(1, 1);
  deltaDD ~ normal(mu_deltaDD, 1 / sqrt(tau_deltaDD));
  omegaDD ~ normal(mu_omegaDD, 1 / sqrt(tau_omegaDD));
  epsilonDD ~ normal(mu_epsilonDD, 1 / sqrt(tau_epsilonDD));
  mu_gammaTR ~ normal(1, 0.5);
  tau_gammaTR ~ gamma(1, 1);
  mu_tauTR ~ normal(1, 0.5);
  tau_tauTR ~ gamma(1, 1);
  mu_kappaTR ~ normal(1, 0.5);
  tau_kappaTR ~ gamma(1, 1);
  mu_varthetaTR ~ normal(1, 0.5);
  tau_varthetaTR ~ gamma(1, 1);
  mu_epsilonTR ~ normal(0.05, 0.1);
  tau_epsilonTR ~ gamma(1, 1);
  gammaTR ~ normal(mu_gammaTR, 1 / sqrt(tau_gammaTR));
  tauTR ~ normal(mu_tauTR, 1 / sqrt(tau_tauTR));
  kappaTR ~ normal(mu_kappaTR, 1 / sqrt(tau_kappaTR));
  varthetaTR ~ normal(mu_varthetaTR, 1 / sqrt(tau_varthetaTR));
  epsilonTR ~ normal(mu_epsilonTR, 1 / sqrt(tau_epsilonTR));
  mu_gammaUT ~ normal(1, 0.5);
  tau_gammaUT ~ gamma(1, 1);
  mu_tauUT ~ normal(1, 0.5);
  tau_tauUT ~ gamma(1, 1);
  mu_kappaUT ~ normal(1, 0.5);
  tau_kappaUT ~ gamma(1, 1);
  mu_varthetaUT ~ normal(1, 0.5);
  tau_varthetaUT ~ gamma(1, 1);
  mu_etaUT ~ normal(0, 0.5);
  tau_etaUT ~ gamma(1, 1);
  mu_epsilonUT ~ normal(0.05, 0.1);
  tau_epsilonUT ~ gamma(1, 1);
  gammaUT ~ normal(mu_gammaUT, 1 / sqrt(tau_gammaUT));
  tauUT ~ normal(mu_tauUT, 1 / sqrt(tau_tauUT));
  kappaUT ~ normal(mu_kappaUT, 1 / sqrt(tau_kappaUT));
  varthetaUT ~ normal(mu_varthetaUT, 1 / sqrt(tau_varthetaUT));
  etaUT ~ normal(mu_etaUT, 1 / sqrt(tau_etaUT));
  epsilonUT ~ normal(mu_epsilonUT, 1 / sqrt(tau_epsilonUT));
  mu_beta0IT ~ normal(0, 0.5);
  tau_beta0IT ~ gamma(1, 1);
  mu_betaRAIT ~ normal(0, 0.5);
  tau_betaRAIT ~ gamma(1, 1);
  mu_betaRRIT ~ normal(0, 0.5);
  tau_betaRRIT ~ gamma(1, 1);
  mu_betaTAIT ~ normal(0, 0.5);
  tau_betaTAIT ~ gamma(1, 1);
  mu_betaTRIT ~ normal(0, 0.5);
  tau_betaTRIT ~ gamma(1, 1);
  mu_epsilonIT ~ normal(0.05, 0.1);
  tau_epsilonIT ~ gamma(1, 1);
  beta0IT ~ normal(mu_beta0IT, 1 / sqrt(tau_beta0IT));
  betaRAIT ~ normal(mu_betaRAIT, 1 / sqrt(tau_betaRAIT));
  betaRRIT ~ normal(mu_betaRRIT, 1 / sqrt(tau_betaRRIT));
  betaTAIT ~ normal(mu_betaTAIT, 1 / sqrt(tau_betaTAIT));
  betaTRIT ~ normal(mu_betaTRIT, 1 / sqrt(tau_betaTRIT));
  epsilonIT ~ normal(mu_epsilonIT, 1 / sqrt(tau_epsilonIT));
  mu_alphaLL ~ normal(0.95, 0.05);
  tau_alphaLL ~ gamma(1, 1);
  mu_alphaSS ~ normal(0.95, 0.05);
  tau_alphaSS ~ gamma(1, 1);
  alphaLL ~ normal(mu_alphaLL, 1 / sqrt(tau_alphaLL));
  alphaSS ~ normal(mu_alphaSS, 1 / sqrt(tau_alphaSS));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      vector[nModels] lp;
      for (m in 1:nModels)
        lp[m] = log_uniform + bernoulli_lpmf(decision[i, j] | theta[i, j, m]);
      target += log_sum_exp(lp);
    }
  }
}
generated quantities {
  array[nParticipants] int z;
  for (i in 1:nParticipants) {
    vector[nModels] score;
    score = rep_vector(0, nModels);
    for (j in 1:nTrials)
      for (m in 1:nModels)
        score[m] += bernoulli_lpmf(decision[i, j] | theta[i, j, m]);
    z[i] = which_max(score);
  }
}
