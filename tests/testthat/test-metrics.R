test_that("perfectly separated scores give AUC 1 and TSS 1", {
  m <- calculate_classifier_metrics(scores_pos = c(0.9, 0.8, 0.95), scores_neg = c(0.1, 0.2, 0.05))
  expect_equal(m$auc_roc, 1)
  expect_equal(m$tss, 1)
  expect_equal(m$ba, 1)
})

test_that("perfectly inverted scores give AUC 0", {
  m <- calculate_classifier_metrics(scores_pos = c(0.1, 0.2), scores_neg = c(0.8, 0.9))
  expect_equal(m$auc_roc, 0)
})

test_that("constant scores do not spuriously report a perfect TSS", {
  # Ties must be grouped: every point sharing a score has to fall on the same side
  # of any threshold. A naive cumulative sweep counts the presences first and
  # reports TSS = 1 on input that carries no signal at all.
  m <- calculate_classifier_metrics(scores_pos = rep(0.5, 20), scores_neg = rep(0.5, 20))
  expect_equal(m$auc_roc, 0.5)
  expect_equal(m$tss, 0)
  expect_equal(m$ba, 0.5)
})

test_that("AUC-ROC matches the Mann-Whitney statistic on a mixed case", {
  pos <- c(3, 1, 4)
  neg <- c(2, 0)
  # pairs where pos > neg: (3>2,3>0), (1>0), (4>2,4>0) = 5 of 6
  expect_equal(calculate_classifier_metrics(pos, neg)$auc_roc, 5 / 6)
})

test_that("empty classes fall back to neutral metrics rather than erroring", {
  m <- calculate_classifier_metrics(scores_pos = numeric(0), scores_neg = c(0.1, 0.2))
  expect_equal(m$auc_roc, 0.5)
  expect_equal(m$tss, 0)
})

test_that("NA scores are dropped before metrics are computed", {
  with_na    <- calculate_classifier_metrics(c(0.9, NA, 0.8), c(0.1, 0.2, NA))
  without_na <- calculate_classifier_metrics(c(0.9, 0.8), c(0.1, 0.2))
  expect_equal(with_na$auc_roc, without_na$auc_roc)
})

test_that("CBI is positive when presences concentrate at high scores", {
  set.seed(1)
  pos <- runif(200, 0.6, 1.0)
  all <- c(pos, runif(800, 0.0, 1.0))
  expect_gt(calculate_cbi(pos, all), 0)
})

test_that("CBI returns 0 when the score range collapses", {
  expect_equal(calculate_cbi(rep(0.5, 10), rep(0.5, 50)), 0)
})
