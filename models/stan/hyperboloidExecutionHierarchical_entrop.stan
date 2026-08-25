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
  real<lower=1e-4, upper=100> mu_kappa;
  real<lower=0.05, upper=25> tau_kappa;
  real<lower=1e-4, upper=100> mu_tau;
  real<lower=0.05, upper=25> tau_tau;
  real<lower=0, upper=20> mu_w;
  real<lower=0.05, upper=25> tau_w;
  vector<lower=1e-4, upper=100>[nParticipants] kappa;
  vector<lower=1e-4, upper=100>[nParticipants] tau;
  vector<lower=0, upper=20>[nParticipants] w;
}

model {
  mu_kappa ~ normal(1, inv_sqrt(4 * muPrecScale));
  tau_kappa ~ gamma(1, 1);
  mu_tau ~ normal(1, inv_sqrt(4 * muPrecScale));
  tau_tau ~ gamma(1, 1);
  mu_w ~ normal(2, inv_sqrt(1 * muPrecScale));
  tau_w ~ gamma(1, 1);
  kappa ~ normal(mu_kappa, inv_sqrt(tau_kappa));
  tau ~ normal(mu_tau, inv_sqrt(tau_tau));
  w ~ normal(mu_w, inv_sqrt(tau_w));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      real valueLL = rLL[i, j] / pow(1 + kappa[i] * tLL[i, j], tau[i]);
      real valueSS = rSS[i, j] / pow(1 + kappa[i] * tSS[i, j], tau[i]);
      decision[i, j] ~ bernoulli_logit(entrop_logit(valueLL - valueSS, w[i]));
    }
  }
}
