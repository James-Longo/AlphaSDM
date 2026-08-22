skip_if_not_installed("dagitty")

# X -> Y confounded by A and B; M mediates X -> Y.
demo_dag <- 'dag {
  X [exposure] Y [outcome]
  A -> X ; A -> Y ; B -> X ; B -> Y ; X -> M -> Y ; X -> Y
}'

sim <- function(n = 2000, seed = 1) {
  set.seed(seed)
  A <- rnorm(n); B <- rnorm(n)
  X <- A + B + rnorm(n)
  M <- 0.5 * X + rnorm(n)
  Y <- 2 * X + M + A + B + rnorm(n)
  data.frame(X, Y, A, B, M)
}

test_that("total-effect adjustment set excludes the mediator", {
  cm <- causal_model(demo_dag, sim())
  adj <- causal_adjustment_sets(cm, "X", "Y")
  expect_true(all(!grepl("M", adj$set)))
  expect_true(any(grepl("A", adj$set) & grepl("B", adj$set)))
})

test_that("sets are ranked by reliability and latent members marked unusable", {
  d <- sim(); d$B <- NULL     # B unmeasured -> the only valid set has a latent node
  cm <- causal_model(demo_dag, d,
                     measurements = list(A = list(reliability = 0.4, source = "idx")))
  adj <- causal_adjustment_sets(cm, "X", "Y")
  expect_false(any(adj$usable))
  expect_error(causal_effect(cm, "X", "Y"), "unmeasured")
})

test_that("causal_effect recovers a known total effect with glm", {
  cm <- causal_model(demo_dag, sim(),
                     measurements = list(A = list(reliability = 0.9, source = "idx"),
                                         B = list(reliability = 0.8, source = "idx")))
  eff <- causal_effect(cm, "X", "Y")
  # total effect = direct 2 + mediated 0.5 * 1 = 2.5
  expect_equal(unname(eff$estimate), 2.5, tolerance = 0.1)
  expect_setequal(eff$adjustment_set, c("A", "B"))
  expect_match(eff$attenuation_note, "0.80")
})

test_that("causal_dag_test flags a broken independence", {
  d <- sim()
  # DAG omitting A -> Y implies A _||_ Y | X, ... which the data violate
  wrong <- 'dag { X [exposure] Y [outcome] A -> X ; B -> X ; B -> Y ; X -> Y }'
  cm <- causal_model(wrong, d[, c("X", "Y", "A", "B")])
  lt <- causal_dag_test(cm)
  expect_true("contradicts_dag" %in% names(lt))
  expect_true(any(lt$contradicts_dag))
})

test_that("causal_report writes a self-contained html", {
  skip_if_not_installed("jsonlite")
  cm <- causal_model(demo_dag, sim())
  f <- withr::local_tempfile(fileext = ".html")
  causal_report(cm, file = f, title = "demo")
  html <- paste(readLines(f), collapse = "\n")
  expect_true(all(vapply(c("X", "Y", "M", "svg", "adjustment_sets"),
                         grepl, logical(1), x = html)))
})
