#include "entrop.stan"

data {
  int<lower=1> nParticipants;
  int<lower=1> nTrials;
  array[nParticipants, nTrials] int<lower=0, upper=1> decision;
  matrix[nParticipants, nTrials] rLL;
  matrix[nParticipants, nTrials] rSS;
  matrix[nParticipants, nTrials] tLL;
  matrix[nParticipants, nTrials] tSS;
  real<lower=0> muPrecScale;
}

parameters {
  real<lower=-2, upper=2> mu_beta0;
  real<lower=0.05, upper=25> tau_beta0;
  real<lower=0, upper=2> mu_betaRA;
  real<lower=0.05, upper=25> tau_betaRA;
  real<lower=0, upper=2> mu_betaRR;
  real<lower=0.05, upper=25> tau_betaRR;
  real<lower=-2, upper=0> mu_betaTA;
  real<lower=0.05, upper=25> tau_betaTA;
  real<lower=-2, upper=0> mu_betaTR;
  real<lower=0.05, upper=25> tau_betaTR;
  real<lower=0, upper=20> mu_w;
  real<lower=0.05, upper=25> tau_w;
  vector<lower=-2, upper=2>[nParticipants] beta0;
  vector<lower=0, upper=2>[nParticipants] betaRA;
  vector<lower=0, upper=2>[nParticipants] betaRR;
  vector<lower=-2, upper=0>[nParticipants] betaTA;
  vector<lower=-2, upper=0>[nParticipants] betaTR;
  vector<lower=0, upper=20>[nParticipants] w;
}

model {
  mu_beta0 ~ normal(0, inv_sqrt(4 * muPrecScale));
  tau_beta0 ~ gamma(1, 1);
  mu_betaRA ~ normal(0, inv_sqrt(4 * muPrecScale));
  tau_betaRA ~ gamma(1, 1);
  mu_betaRR ~ normal(0, inv_sqrt(4 * muPrecScale));
  tau_betaRR ~ gamma(1, 1);
  mu_betaTA ~ normal(0, inv_sqrt(4 * muPrecScale));
  tau_betaTA ~ gamma(1, 1);
  mu_betaTR ~ normal(0, inv_sqrt(4 * muPrecScale));
  tau_betaTR ~ gamma(1, 1);
  mu_w ~ normal(2, inv_sqrt(1 * muPrecScale));
  tau_w ~ gamma(1, 1);
  beta0 ~ normal(mu_beta0, inv_sqrt(tau_beta0));
  betaRA ~ normal(mu_betaRA, inv_sqrt(tau_betaRA));
  betaRR ~ normal(mu_betaRR, inv_sqrt(tau_betaRR));
  betaTA ~ normal(mu_betaTA, inv_sqrt(tau_betaTA));
  betaTR ~ normal(mu_betaTR, inv_sqrt(tau_betaTR));
  w ~ normal(mu_w, inv_sqrt(tau_w));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      real dRewardA = betaRA[i] * (rLL[i, j] - rSS[i, j]);
      real dRewardR = betaRR[i] * ((rLL[i, j] - rSS[i, j]) / (0.5 * (rLL[i, j] + rSS[i, j])));
      real dTimeA = betaTA[i] * (tLL[i, j] - tSS[i, j]);
      real dTimeR = betaTR[i] * ((tLL[i, j] - tSS[i, j]) / (0.5 * (tLL[i, j] + tSS[i, j])));
      decision[i, j] ~ bernoulli_logit(
        entrop_logit(beta0[i] + dRewardA + dRewardR + dTimeA + dTimeR, w[i]));
    }
  }
}
