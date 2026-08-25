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
  real<lower=-2000, upper=2000> mu_delta;
  real<lower=0.05, upper=25> tau_delta;
  real<lower=1e-6, upper=1 - 1e-6> mu_omega;
  real<lower=0.05, upper=25> tau_omega;
  real<lower=0, upper=20> mu_w;
  real<lower=0.05, upper=25> tau_w;
  vector<lower=-2000, upper=2000>[nParticipants] delta;
  vector<lower=1e-6, upper=1 - 1e-6>[nParticipants] omega;
  vector<lower=0, upper=20>[nParticipants] w;
}

model {
  mu_delta ~ normal(0, inv_sqrt(1e-6 * muPrecScale));
  tau_delta ~ gamma(1, 1);
  mu_omega ~ normal(0.5, inv_sqrt(4 * muPrecScale));
  tau_omega ~ gamma(1, 1);
  mu_w ~ normal(2, inv_sqrt(1 * muPrecScale));
  tau_w ~ gamma(1, 1);
  delta ~ normal(mu_delta, inv_sqrt(tau_delta));
  omega ~ normal(mu_omega, inv_sqrt(tau_omega));
  w ~ normal(mu_w, inv_sqrt(tau_w));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      real dRewardR = omega[i] * (rLL[i, j] - rSS[i, j]);
      real dTimeR = (1 - omega[i]) * (tLL[i, j] - tSS[i, j]);
      decision[i, j] ~ bernoulli_logit(entrop_logit(dRewardR - dTimeR - delta[i], w[i]));
    }
  }
}
