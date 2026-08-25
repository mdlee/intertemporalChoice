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
  real<lower=1e-4, upper=100> mu_kappa;
  real<lower=0.05, upper=25> tau_kappa;
  real<lower=1e-4, upper=100> mu_tau;
  real<lower=0.05, upper=25> tau_tau;
  real<lower=0, upper=0.5> mu_sigma;
  real<lower=0.05, upper=25> tau_sigma;
  vector<lower=1e-4, upper=100>[nParticipants] kappa;
  vector<lower=1e-4, upper=100>[nParticipants] tau;
  vector<lower=1e-4, upper=0.499>[nParticipants] sigma;
}
model {
  mu_kappa ~ normal(1, 0.5);
  tau_kappa ~ gamma(1, 1);
  mu_tau ~ normal(1, 0.5);
  tau_tau ~ gamma(1, 1);
  mu_sigma ~ normal(0.05, 0.1);
  tau_sigma ~ gamma(1, 1);
  kappa ~ normal(mu_kappa, 1 / sqrt(tau_kappa));
  tau ~ normal(mu_tau, 1 / sqrt(tau_tau));
  sigma ~ normal(mu_sigma, 1 / sqrt(tau_sigma));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      real vLL;
      real vSS;
      vLL = rLL[i, j] / pow(1 + kappa[i] * tLL[i, j], tau[i]);
      vSS = rSS[i, j] / pow(1 + kappa[i] * tSS[i, j], tau[i]);
      decision[i, j] ~ bernoulli(probit_theta(vLL - vSS, sigma[i]));
    }
  }
}
