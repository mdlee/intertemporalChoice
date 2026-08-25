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
  real<lower=1e-4, upper=100> mu_gamma;
  real<lower=0.05, upper=25> tau_gamma;
  real<lower=1e-4, upper=100> mu_tau;
  real<lower=0.05, upper=25> tau_tau;
  real<lower=1e-4, upper=100> mu_kappa;
  real<lower=0.05, upper=25> tau_kappa;
  real<lower=1e-4, upper=100> mu_vartheta;
  real<lower=0.05, upper=25> tau_vartheta;
  real<lower=0, upper=0.5> mu_epsilon;
  real<lower=0.05, upper=25> tau_epsilon;
  vector<lower=1e-4, upper=100>[nParticipants] gamma;
  vector<lower=1e-4, upper=100>[nParticipants] tau;
  vector<lower=1e-4, upper=100>[nParticipants] kappa;
  vector<lower=1e-4, upper=100>[nParticipants] vartheta;
  vector<lower=1e-4, upper=0.499>[nParticipants] epsilon;
}
model {
  mu_gamma ~ normal(1, 0.5);
  tau_gamma ~ gamma(1, 1);
  mu_tau ~ normal(1, 0.5);
  tau_tau ~ gamma(1, 1);
  mu_kappa ~ normal(1, 0.5);
  tau_kappa ~ gamma(1, 1);
  mu_vartheta ~ normal(1, 0.5);
  tau_vartheta ~ gamma(1, 1);
  mu_epsilon ~ normal(0.05, 0.1);
  tau_epsilon ~ gamma(1, 1);
  gamma ~ normal(mu_gamma, 1 / sqrt(tau_gamma));
  tau ~ normal(mu_tau, 1 / sqrt(tau_tau));
  kappa ~ normal(mu_kappa, 1 / sqrt(tau_kappa));
  vartheta ~ normal(mu_vartheta, 1 / sqrt(tau_vartheta));
  epsilon ~ normal(mu_epsilon, 1 / sqrt(tau_epsilon));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      real vLL;
      real vSS;
      real wLL;
      real wSS;
      real qV;
      real qT;
      vLL = (1 / gamma[i]) * log(1 + gamma[i] * rLL[i, j]);
      vSS = (1 / gamma[i]) * log(1 + gamma[i] * rSS[i, j]);
      wLL = (1 / tau[i]) * log(1 + tau[i] * tLL[i, j]);
      wSS = (1 / tau[i]) * log(1 + tau[i] * tSS[i, j]);
      qV = vLL - vSS;
      qT = tradeoff_qtime(wLL, wSS, kappa[i], vartheta[i]);
      decision[i, j] ~ bernoulli(execution_theta(qV - qT, epsilon[i]));
    }
  }
}
