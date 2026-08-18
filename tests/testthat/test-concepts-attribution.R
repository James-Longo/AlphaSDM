test_that("selection_scores recovers a planted enrichment", {
  set.seed(2)
  m <- 6L
  bands <- concept_band_names(m)
  bg <- matrix(stats::rexp(500 * m), 500, m, dimnames = list(NULL, bands))
  pres <- matrix(stats::rexp(300 * m), 300, m, dimnames = list(NULL, bands))
  pres[, 3] <- pres[, 3] * 4        # the planted concept

  s <- selection_scores(pres, bg, boot_reps = 200L)
  expect_equal(s$concept_id, bands)
  expect_equal(which.max(s$selection_ratio), 3L)
  expect_gt(s$selection_lo[3], 2)   # CI excludes no-enrichment
  expect_true(all(s$selection_lo <= s$selection_ratio + 1e-8))
  expect_true(all(s$selection_hi >= s$selection_ratio - 1e-8))
})

test_that("a single supporting pixel cannot hold an infinite lower bound", {
  # One presence pixel, never in background: the ratio is legitimately Inf,
  # but resamples that miss that pixel are 0/0 = zero evidence, so the
  # bootstrap lower bound must collapse rather than stay Inf.
  set.seed(4)
  m <- 2L
  bands <- concept_band_names(m)
  pres <- matrix(0, 100, m, dimnames = list(NULL, bands))
  bg <- matrix(0, 100, m, dimnames = list(NULL, bands))
  pres[, 1] <- stats::rexp(100); bg[, 1] <- stats::rexp(100)
  pres[1, 2] <- 0.5                       # the one-pixel concept
  s <- selection_scores(pres, bg, boot_reps = 400L)
  expect_true(is.infinite(s$selection_ratio[2]))
  expect_equal(s$selection_lo[2], 0)
  expect_equal(s$pres_active[2], 0.01)
  expect_equal(s$bg_active[2], 0)
})

test_that("selection_scores reports silent concepts honestly", {
  m <- 3L
  bands <- concept_band_names(m)
  bg <- matrix(1, 50, m, dimnames = list(NULL, bands))
  pres <- matrix(1, 50, m, dimnames = list(NULL, bands))
  bg[, 2] <- 0                       # active at presences only -> Inf
  bg[, 3] <- 0; pres[, 3] <- 0       # silent everywhere -> NaN
  s <- selection_scores(pres, bg, boot_reps = 10L)
  expect_equal(s$selection_ratio[1], 1)
  expect_true(is.infinite(s$selection_ratio[2]))
  expect_true(is.nan(s$selection_ratio[3]))
})

test_that("selection_scores drops rows with missing activations", {
  m <- 3L
  bands <- concept_band_names(m)
  pres <- matrix(1, 10, m, dimnames = list(NULL, bands))
  bg <- matrix(1, 10, m, dimnames = list(NULL, bands))
  pres[1, 2] <- NA
  s <- selection_scores(pres, bg, boot_reps = 10L)
  expect_false(anyNA(s$selection_ratio))
  expect_error(selection_scores(pres[0, , drop = FALSE], bg, boot_reps = 10L),
               "No unmasked")
})

test_that("projection_scores recovers a planted direction", {
  sae <- make_test_sae(m = 32L, k = 8L)
  j <- 7L
  beta <- sae$W_dec[j, ] * 2.5       # linear model aligned with concept j
  p <- projection_scores(beta, sae)
  expect_equal(which.max(p$projection_cos), j)
  expect_equal(p$projection_cos[j], 1)
  # Random unit rows in 64-d are near-orthogonal, so nothing else comes close.
  expect_true(all(p$projection_cos[-j] < 0.6))
  expect_error(projection_scores(beta[-1], sae), "64-d")
})

test_that("ablating the planted concept collapses a one-concept linear score", {
  # A synthetic SDM that scores by dot product with one decoder direction:
  # removing that concept must erase the signal, removing another must not.
  # Tied weights, so encoding a decoder-built input recovers the planted
  # activation (near-orthogonal random unit rows make the Gram near-identity).
  rnd <- make_test_sae(m = 16L, k = 4L)
  sae <- new_sae(W_enc = rnd$W_dec, b_enc = rep(0, 16), W_dec = rnd$W_dec,
                 b_pre = rep(0, 64), k = 4L)
  set.seed(3)
  j <- 5L
  n <- 50L
  acts <- matrix(0, n, sae$m)
  acts[, j] <- stats::rexp(n) + 1
  X <- sae_reconstruct(acts, sae)
  w <- sae$W_dec[j, ]
  score <- function(M) as.numeric(M %*% w)

  base <- score(X)
  drop_j <- mean(base - score(sae_ablate(X, sae, j)))
  others <- setdiff(which(colSums(sae_encode(X, sae)) > 0), j)
  drops_other <- vapply(others, function(o)
    mean(base - score(sae_ablate(X, sae, o))), numeric(1))
  expect_gt(drop_j, 0)
  expect_true(all(drop_j > drops_other))
})

test_that("derive_concepts validates its sdm argument", {
  expect_error(derive_concepts(list(methods = "rf")), "evaluate_models")
  expect_error(derive_concepts(list(models = list(), context = NULL)), "evaluate_models")
})

test_that("the print method flags unlabeled concepts", {
  df <- data.frame(concept_id = c("C001", "C002"),
                   label = c("tree cover, high canopy", NA),
                   coherence = c(0.8, NA),
                   selection_ratio = c(3.1, 2.2),
                   stringsAsFactors = FALSE)
  class(df) <- c("alphasdm_concepts", "data.frame")
  attr(df, "methods") <- "selection"
  attr(df, "model") <- "rf"
  attr(df, "top_n") <- 2L
  out <- paste(utils::capture.output(print(df)), collapse = "\n")
  expect_match(out, "unlabeled latent")
  expect_match(out, "do not over-interpret")
})

test_that("example location strings round-trip", {
  df <- parse_example_locations("-73.60000,44.50000;10.10000,-2.20000")
  expect_equal(nrow(df), 2L)
  expect_equal(df$longitude, c(-73.6, 10.1))
  expect_equal(df$latitude, c(44.5, -2.2))
  expect_null(parse_example_locations(""))
  expect_null(parse_example_locations(NA))
})
