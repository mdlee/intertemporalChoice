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
  real<lower=1e-4, upper=100> mu_gamma;
  real<lower=0.05, upper=25> tau_gamma;
  real<lower=1e-4, upper=100> mu_tau;
  real<lower=0.05, upper=25> tau_tau;
  real<lower=1e-4, upper=100> mu_kappa;
  real<lower=0.05, upper=25> tau_kappa;
  real<lower=1e-4, upper=100> mu_vartheta;
  real<lower=0.05, upper=25> tau_vartheta;
  real<lower=0, upper=500> mu_eta;
  real<lower=0.05, upper=25> tau_eta;
  real<lower=0, upper=20> mu_w;
  real<lower=0.05, upper=25> tau_w;
  vector<lower=1e-4, upper=100>[nParticipants] gamma;
  vector<lower=1e-4, upper=100>[nParticipants] tau;
  vector<lower=1e-4, upper=100>[nParticipants] kappa;
  vector<lower=1e-4, upper=100>[nParticipants] vartheta;
  vector<lower=0, upper=500>[nParticipants] eta;
  vector<lower=0, upper=20>[nParticipants] w;
}

model {
  mu_gamma ~ normal(1, inv_sqrt(4 * muPrecScale));
  tau_gamma ~ gamma(1, 1);
  mu_tau ~ normal(1, inv_sqrt(4 * muPrecScale));
  tau_tau ~ gamma(1, 1);
  mu_kappa ~ normal(1, inv_sqrt(4 * muPrecScale));
  tau_kappa ~ gamma(1, 1);
  mu_vartheta ~ normal(1, inv_sqrt(4 * muPrecScale));
  tau_vartheta ~ gamma(1, 1);
  mu_eta ~ normal(0, inv_sqrt(4 * muPrecScale));
  tau_eta ~ gamma(1, 1);
  mu_w ~ normal(2, inv_sqrt(1 * muPrecScale));
  tau_w ~ gamma(1, 1);
  gamma ~ normal(mu_gamma, inv_sqrt(tau_gamma));
  tau ~ normal(mu_tau, inv_sqrt(tau_tau));
  kappa ~ normal(mu_kappa, inv_sqrt(tau_kappa));
  vartheta ~ normal(mu_vartheta, inv_sqrt(tau_vartheta));
  eta ~ normal(mu_eta, inv_sqrt(tau_eta));
  w ~ normal(mu_w, inv_sqrt(tau_w));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      real rLLnet = fmax(rLL[i, j] - eta[i], 1e-6);
      real valueLL = (1 / gamma[i]) * log(1 + gamma[i] * rLLnet);
      real valueSS = (1 / gamma[i]) * log(1 + gamma[i] * rSS[i, j]);
      real weightLL = (1 / tau[i]) * log(1 + tau[i] * tLL[i, j]);
      real weightSS = (1 / tau[i]) * log(1 + tau[i] * tSS[i, j]);
      real qValue = valueLL - valueSS;
      real dW = fmax(weightLL - weightSS, 0);
      real qTime = kappa[i] * log(1 + pow(dW / vartheta[i], vartheta[i]));
      decision[i, j] ~ bernoulli_logit(entrop_logit(qValue - qTime, w[i]));
    }
  }
}
