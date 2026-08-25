#include "common_functions.stan"

data {
  int<lower=1> nParticipants;
  int<lower=1> nTrials;
  matrix[nParticipants, nTrials] rLL;
  matrix[nParticipants, nTrials] rSS;
  matrix[nParticipants, nTrials] tLL;
  matrix[nParticipants, nTrials] tSS;
  array[nParticipants, nTrials] int<lower=0, upper=1> decision;
}
parameters {
  real<lower=-2, upper=2> mu_delta;
  real<lower=0.05, upper=25> tau_delta;
  real<lower=0, upper=0.5> mu_sigma;
  real<lower=0.05, upper=25> tau_sigma;
  vector<lower=-2, upper=2>[nParticipants] delta;
  vector<lower=1e-4, upper=0.499>[nParticipants] sigma;
}
model {
  mu_delta ~ normal(0, 0.5);
  tau_delta ~ gamma(1, 1);
  mu_sigma ~ normal(0.05, 0.1);
  tau_sigma ~ gamma(1, 1);
  delta ~ normal(mu_delta, 1 / sqrt(tau_delta));
  sigma ~ normal(mu_sigma, 1 / sqrt(tau_sigma));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      real dr;
      real dt;
      dr = (rLL[i, j] - rSS[i, j]) / rLL[i, j];
      dt = (tLL[i, j] - tSS[i, j]) / tLL[i, j];
      decision[i, j] ~ bernoulli(probit_theta(dr - dt - delta[i], sigma[i]));
    }
  }
}
