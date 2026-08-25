functions {
// Shared decision rules (match JAGS execution / probit link)
real execution_theta(real diff, real eps) {
  return diff > 0.0 ? (1.0 - eps) : eps;
}
real probit_theta(real latent, real sigma) {
  return Phi(latent / sigma);
}
real tradeoff_qtime(real wLL, real wSS, real kappa, real vartheta) {
  real dw;
  dw = wLL - wSS;
  return kappa * log(1.0 + pow(fmax(dw / vartheta, 0.0), vartheta));
}
real mix_ex_execution(real rLL, real rSS, real tLL, real tSS, real kappa, real eps) {
  return execution_theta(rLL * exp(-kappa * tLL) - rSS * exp(-kappa * tSS), eps);
}
real mix_hc_execution(real rLL, real rSS, real tLL, real tSS, real kappa, real eps) {
  return execution_theta(rLL / (1 + kappa * tLL) - rSS / (1 + kappa * tSS), eps);
}
real mix_hy_execution(real rLL, real rSS, real tLL, real tSS, real kappa, real tau, real eps) {
  return execution_theta(
    rLL / pow(1 + kappa * tLL, tau) - rSS / pow(1 + kappa * tSS, tau), eps);
}
real mix_pd_execution(real rLL, real rSS, real tLL, real tSS, real delta, real eps) {
  return execution_theta(
    (rLL - rSS) / rLL - (tLL - tSS) / tLL - delta, eps);
}
real mix_dd_execution(real rLL, real rSS, real tLL, real tSS, real delta, real omega, real eps) {
  return execution_theta(
    omega * (rLL - rSS) - (1 - omega) * (tLL - tSS) - delta, eps);
}
real mix_tr_execution(real rLL, real rSS, real tLL, real tSS,
    real gamma, real tau, real kappa, real vartheta, real eps) {
  real vLL;
  real vSS;
  real wLL;
  real wSS;
  vLL = (1 / gamma) * log(1 + gamma * rLL);
  vSS = (1 / gamma) * log(1 + gamma * rSS);
  wLL = (1 / tau) * log(1 + tau * tLL);
  wSS = (1 / tau) * log(1 + tau * tSS);
  return execution_theta(vLL - vSS - tradeoff_qtime(wLL, wSS, kappa, vartheta), eps);
}
real mix_ut_execution(real rLL, real rSS, real tLL, real tSS,
    real gamma, real tau, real kappa, real vartheta, real eta, real eps) {
  real vLL;
  real vSS;
  real wLL;
  real wSS;
  vLL = (1 / gamma) * log(1 + gamma * (rLL - eta));
  vSS = (1 / gamma) * log(1 + gamma * rSS);
  wLL = (1 / tau) * log(1 + tau * tLL);
  wSS = (1 / tau) * log(1 + tau * tSS);
  return execution_theta(vLL - vSS - tradeoff_qtime(wLL, wSS, kappa, vartheta), eps);
}
real mix_it_execution(real rLL, real rSS, real tLL, real tSS,
    real b0, real bRA, real bRR, real bTA, real bTR, real eps) {
  real score;
  score = b0 + bRA * (rLL - rSS)
    + bRR * ((rLL - rSS) / (0.5 * (rLL + rSS)))
    + bTA * (tLL - tSS)
    + bTR * ((tLL - tSS) / (0.5 * (tLL + tSS)));
  return execution_theta(score, eps);
}
real mix_ex_probit(real rLL, real rSS, real tLL, real tSS, real kappa, real sigma) {
  return probit_theta(rLL * exp(-kappa * tLL) - rSS * exp(-kappa * tSS), sigma);
}
real mix_hc_probit(real rLL, real rSS, real tLL, real tSS, real kappa, real sigma) {
  return probit_theta(rLL / (1 + kappa * tLL) - rSS / (1 + kappa * tSS), sigma);
}
real mix_hy_probit(real rLL, real rSS, real tLL, real tSS, real kappa, real tau, real sigma) {
  return probit_theta(
    rLL / pow(1 + kappa * tLL, tau) - rSS / pow(1 + kappa * tSS, tau), sigma);
}
real mix_pd_probit(real rLL, real rSS, real tLL, real tSS, real delta, real sigma) {
  return probit_theta(
    (rLL - rSS) / rLL - (tLL - tSS) / tLL - delta, sigma);
}
real mix_dd_probit(real rLL, real rSS, real tLL, real tSS, real delta, real omega, real sigma) {
  return probit_theta(
    omega * (rLL - rSS) - (1 - omega) * (tLL - tSS) - delta, sigma);
}
real mix_tr_probit(real rLL, real rSS, real tLL, real tSS,
    real gamma, real tau, real kappa, real vartheta, real sigma) {
  real vLL;
  real vSS;
  real wLL;
  real wSS;
  vLL = (1 / gamma) * log(1 + gamma * rLL);
  vSS = (1 / gamma) * log(1 + gamma * rSS);
  wLL = (1 / tau) * log(1 + tau * tLL);
  wSS = (1 / tau) * log(1 + tau * tSS);
  return probit_theta(vLL - vSS - tradeoff_qtime(wLL, wSS, kappa, vartheta), sigma);
}
real mix_ut_probit(real rLL, real rSS, real tLL, real tSS,
    real gamma, real tau, real kappa, real vartheta, real eta, real sigma) {
  real vLL;
  real vSS;
  real wLL;
  real wSS;
  vLL = (1 / gamma) * log(1 + gamma * (rLL - eta));
  vSS = (1 / gamma) * log(1 + gamma * rSS);
  wLL = (1 / tau) * log(1 + tau * tLL);
  wSS = (1 / tau) * log(1 + tau * tSS);
  return probit_theta(vLL - vSS - tradeoff_qtime(wLL, wSS, kappa, vartheta), sigma);
}
real mix_it_probit(real rLL, real rSS, real tLL, real tSS,
    real b0, real bRA, real bRR, real bTA, real bTR, real sigma) {
  real score;
  score = b0 + bRA * (rLL - rSS)
    + bRR * ((rLL - rSS) / (0.5 * (rLL + rSS)))
    + bTA * (tLL - tSS)
    + bTR * ((tLL - tSS) / (0.5 * (tLL + tSS)));
  return probit_theta(score, sigma);
}
}
