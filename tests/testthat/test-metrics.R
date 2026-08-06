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

test_that("a request Earth Engine refuses is recognised as one to escalate", {
  # The classifier in export_image_tiled decides whether to move a tile into the
  # batch system or give up. It missed "Computation timed out" once, which is the
  # single most common way a fine-scale tile fails, so the cases are pinned here.
  escalate <- function(msg) {
    is_gee_timeout(msg) || grepl("400", msg, fixed = TRUE) ||
      grepl("Timeout of", msg, fixed = TRUE)
  }
  expect_true(escalate("ee.ee_exception.EEException: Computation timed out."))
  expect_true(escalate("User memory limit exceeded."))
  expect_true(escalate("cannot open URL: HTTP status was '400 Bad Request'"))
  expect_true(escalate("Timeout of 900 seconds was reached"))
  expect_true(escalate("Too many concurrent aggregations"))
  # A network problem is not fixed by computing the image somewhere else.
  expect_false(escalate("cannot open the connection"))
  expect_false(escalate("Could not resolve host: earthengine.googleapis.com"))
})
