functions {
  // Grünwald entropification of a hard integration rule:
  // P(choose LL | rule favors LL) = inv_logit(w)
  // P(choose LL | rule favors SS) = inv_logit(-w)
  real entrop_logit(real advantage, real w) {
    return w * (2 * step(advantage) - 1);
  }
}
