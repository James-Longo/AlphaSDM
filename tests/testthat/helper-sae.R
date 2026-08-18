# Synthetic SAE weights for the concept tests. Small m so the tests are fast;
# unit-norm decoder rows, matching what training ships.
make_test_sae <- function(m = 16L, k = 4L, seed = 1L) {
  set.seed(seed)
  W_dec <- matrix(stats::rnorm(m * 64), m, 64)
  W_dec <- W_dec / sqrt(rowSums(W_dec^2))
  new_sae(W_enc = matrix(stats::rnorm(m * 64, sd = 0.3), m, 64),
          b_enc = stats::rnorm(m, sd = 0.1),
          W_dec = W_dec,
          b_pre = stats::rnorm(64, sd = 0.05),
          k = k)
}
