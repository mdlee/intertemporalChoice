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
  real<lower=-2, upper=2> mu_beta0;
  real<lower=0.05, upper=25> tau_beta0;
  real<lower=-2, upper=2> mu_betaRA;
  real<lower=0.05, upper=25> tau_betaRA;
  real<lower=-2, upper=2> mu_betaRR;
  real<lower=0.05, upper=25> tau_betaRR;
  real<lower=-2, upper=2> mu_betaTA;
  real<lower=0.05, upper=25> tau_betaTA;
  real<lower=-2, upper=2> mu_betaTR;
  real<lower=0.05, upper=25> tau_betaTR;
  real<lower=0, upper=0.5> mu_sigma;
  real<lower=0.05, upper=25> tau_sigma;
  vector<lower=-2, upper=2>[nParticipants] beta0;
  vector<lower=-2, upper=2>[nParticipants] betaRA;
  vector<lower=-2, upper=2>[nParticipants] betaRR;
  vector<lower=-2, upper=2>[nParticipants] betaTA;
  vector<lower=-2, upper=2>[nParticipants] betaTR;
  vector<lower=1e-4, upper=0.499>[nParticipants] sigma;
}
model {
  mu_beta0 ~ normal(0, 0.5);
  tau_beta0 ~ gamma(1, 1);
  mu_betaRA ~ normal(0, 0.5);
  tau_betaRA ~ gamma(1, 1);
  mu_betaRR ~ normal(0, 0.5);
  tau_betaRR ~ gamma(1, 1);
  mu_betaTA ~ normal(0, 0.5);
  tau_betaTA ~ gamma(1, 1);
  mu_betaTR ~ normal(0, 0.5);
  tau_betaTR ~ gamma(1, 1);
  mu_sigma ~ normal(0.05, 0.1);
  tau_sigma ~ gamma(1, 1);
  beta0 ~ normal(mu_beta0, 1 / sqrt(tau_beta0));
  betaRA ~ normal(mu_betaRA, 1 / sqrt(tau_betaRA));
  betaRR ~ normal(mu_betaRR, 1 / sqrt(tau_betaRR));
  betaTA ~ normal(mu_betaTA, 1 / sqrt(tau_betaTA));
  betaTR ~ normal(mu_betaTR, 1 / sqrt(tau_betaTR));
  sigma ~ normal(mu_sigma, 1 / sqrt(tau_sigma));
  for (i in 1:nParticipants) {
    for (j in 1:nTrials) {
      real score;
      score = beta0[i]
        + betaRA[i] * (rLL[i, j] - rSS[i, j])
        + betaRR[i] * ((rLL[i, j] - rSS[i, j]) / (0.5 * (rLL[i, j] + rSS[i, j])))
        + betaTA[i] * (tLL[i, j] - tSS[i, j])
        + betaTR[i] * ((tLL[i, j] - tSS[i, j]) / (0.5 * (tLL[i, j] + tSS[i, j])));
      decision[i, j] ~ bernoulli(probit_theta(score, sigma[i]));
    }
  }
}
