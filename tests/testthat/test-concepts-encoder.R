test_that("concept band names are zero-padded and unique", {
  expect_equal(concept_band_names(512)[1], "C001")
  expect_equal(concept_band_names(512)[512], "C512")
  expect_equal(concept_band_names(1024)[1], "C0001")
  expect_equal(anyDuplicated(concept_band_names(2048)), 0L)
})

test_that("validate_sae rejects malformed weights", {
  sae <- make_test_sae()
  expect_error(new_sae(sae$W_enc[, 1:32], sae$b_enc, sae$W_dec, sae$b_pre, sae$k),
               "W_enc")
  expect_error(new_sae(sae$W_enc, sae$b_enc[-1], sae$W_dec, sae$b_pre, sae$k),
               "b_enc")
  expect_error(new_sae(sae$W_enc, sae$b_enc, sae$W_dec, sae$b_pre[-1], sae$k),
               "b_pre")
  expect_error(new_sae(sae$W_enc, sae$b_enc, sae$W_dec, sae$b_pre, 0L), "k")
  bad <- sae$W_enc; bad[1, 1] <- NA
  expect_error(new_sae(bad, sae$b_enc, sae$W_dec, sae$b_pre, sae$k), "finite")
})

test_that("sae_encode applies ReLU and keeps exactly the top k", {
  # Identity-block encoder so activations are directly readable from the input.
  m <- 8L
  W_enc <- cbind(diag(m), matrix(0, m, 64 - m))
  sae <- new_sae(W_enc, rep(0, m), matrix(1 / 8, m, 64), rep(0, 64), k = 2L)
  x <- matrix(0, 1, 64)
  x[1, 1:8] <- c(5, 3, -2, 1, 0.5, 4, -1, 0.1)
  a <- sae_encode(x, sae)
  expect_equal(as.numeric(a), c(5, 0, 0, 0, 0, 4, 0, 0))
})

test_that("TopK keeps ties at the threshold", {
  m <- 8L
  W_enc <- cbind(diag(m), matrix(0, m, 64 - m))
  sae <- new_sae(W_enc, rep(0, m), matrix(1 / 8, m, 64), rep(0, 64), k = 2L)
  x <- matrix(0, 1, 64)
  x[1, 1:8] <- c(5, 3, 3, 1, 0, 0, 0, 0)
  a <- sae_encode(x, sae)
  expect_equal(as.numeric(a), c(5, 3, 3, 0, 0, 0, 0, 0))
})

test_that("fewer than k positive activations survive untouched", {
  m <- 8L
  W_enc <- cbind(diag(m), matrix(0, m, 64 - m))
  sae <- new_sae(W_enc, rep(0, m), matrix(1 / 8, m, 64), rep(0, 64), k = 4L)
  x <- matrix(0, 1, 64)
  x[1, 1] <- 2
  a <- sae_encode(x, sae)
  expect_equal(as.numeric(a), c(2, rep(0, 7)))
})

test_that("encode respects the pre-encoder bias and reconstruction inverts it", {
  sae <- make_test_sae(m = 16L, k = 16L)   # k = m: no TopK, pure linear + ReLU
  X <- matrix(stats::rnorm(5 * 64, sd = 0.3), 5, 64)
  A <- sae_encode(X, sae)
  manual <- sweep(X, 2L, sae$b_pre) %*% t(sae$W_enc)
  manual <- pmax(sweep(manual, 2L, sae$b_enc, `+`), 0)
  expect_equal(unname(A), unname(manual))
  expect_equal(dim(sae_reconstruct(A, sae)), c(5L, 64L))
})

test_that("sae_ablate removes exactly the a_j d_j term", {
  sae <- make_test_sae()
  X <- matrix(stats::rnorm(10 * 64, sd = 0.3), 10, 64)
  A <- sae_encode(X, sae)
  j <- which.max(colSums(A > 0))
  expect_equal(sae_ablate(X, sae, j),
               X - A[, j, drop = FALSE] %*% sae$W_dec[j, , drop = FALSE],
               ignore_attr = TRUE)
  # Pixels where concept j is inactive are untouched.
  off <- A[, j] == 0
  expect_equal(sae_ablate(X, sae, j)[off, ], X[off, ], ignore_attr = TRUE)
})

test_that("validate_sae strips dimnames inherited from training data", {
  # Named vectors become Python dicts through reticulate, which ee.Array
  # rejects; weights trained on named matrices must come out clean.
  sae <- make_test_sae()
  colnames(sae$W_enc) <- sprintf("A%02d", 0:63)
  names(sae$b_pre) <- sprintf("A%02d", 0:63)
  names(sae$b_enc) <- paste0("L", seq_len(sae$m))
  out <- validate_sae(sae)
  expect_null(dimnames(out$W_enc))
  expect_null(names(out$b_pre))
  expect_null(names(out$b_enc))
})

test_that("load_sae_weights errors informatively when nothing is installed", {
  # No weights ship while the concept layer is experimental; a missing file
  # must say so rather than fail on system.file's empty string.
  path <- system.file("extdata", "sae_weights_v1.rds", package = "AlphaSDM")
  skip_if(nzchar(path), "weights are installed")
  expect_error(load_sae_weights(), "dev/02_train_sae.R")
})
