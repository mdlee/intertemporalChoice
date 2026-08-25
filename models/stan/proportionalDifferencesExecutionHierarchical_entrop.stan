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
  real<lower=-2, upper=2> mu_delta;
  real<lower=0.05, upper=25> tau_delta;
  real<lower=0, upper=20> mu_w;
  real<lower=0.05, upper=25> tau_w;
  vector<lower=-2, upper=2>[nParticipants] delta;
  vector<lower=0, upper=20>[nParticipants] w;
}

model {
  mu_delta ~ normal(0, inv_sqrt(4 * muPrecScale));
  tau_delta ~ gamma(1, 1);
  mu_w ~ normal(2, inv_sqrt(1 * muPrecScale));
  tau_w ~ gamma(1, 1);
  delta ~ normal(mu_delta, inv_sqrt(tau_delta));
  w ~ normal(mu_w, inv_sqrt(tau_w));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      real dRewardR = (rLL[i, j] - rSS[i, j]) / rLL[i, j];
      real dTimeR = (tLL[i, j] - tSS[i, j]) / tLL[i, j];
      decision[i, j] ~ bernoulli_logit(entrop_logit(dRewardR - dTimeR - delta[i], w[i]));
    }
  }
}
