test_that("envelope threshold is bias-corrected above the observed max", {
  set.seed(1)
  P <- matrix(rnorm(500 * 64), 500, 64)
  env <- estimate_embedding_envelope(P)
  d <- sqrt(rowSums(sweep(sweep(P, 2, env$mu), 2, env$sd, "/")^2))
  expect_gte(env$threshold, max(d))
  expect_true(all(d <= env$threshold))
})

test_that("envelope correction widens with fewer presences", {
  set.seed(2)
  P <- matrix(rnorm(2000 * 64), 2000, 64)
  gap <- function(n) {
    Pi <- P[seq_len(n), , drop = FALSE]
    env <- estimate_embedding_envelope(Pi, seed = 3L)
    d <- sqrt(rowSums(sweep(sweep(Pi, 2, env$mu), 2, env$sd, "/")^2))
    (env$threshold - max(d)) / max(d)
  }
  expect_gt(gap(30), gap(1000))
})

test_that("legacy path requires no presence data and keeps the FC contract", {
  skip_if_not(identical(Sys.getenv("ALPHASDM_LIVE_TESTS"), "1"),
              "live GEE test")
  ensure_gee_authenticated()
  ee <- reticulate::import("ee")
  reg <- ee$Geometry$Rectangle(c(-66.5, 46.5, -66.4, 46.6))
  bg <- generate_background_fc_gee(2023L, 25L, reg)
  expect_true(is.na(bg$radius_m))
  info <- bg$fc$limit(3L)$getInfo()
  expect_equal(info$features[[1]]$properties$present, 0L)
})

test_that("format_data standardizes presence-only data with a directive", {
  d <- data.frame(longitude = c(-66.5, -66.6), latitude = c(47.2, 47.3),
                  year = c(2022, 2023))
  ok <- format_data(d, coords = c("longitude", "latitude"), year = "year")
  expect_equal(nrow(ok), 2L)
  expect_true(all(ok$present == 1))
})

test_that("modelling rejects presence-only data with a directive", {
  d <- data.frame(longitude = c(-66.5, -66.6), latitude = c(47.2, 47.3),
                  year = c(2022, 2023), present = 1)
  expect_error(evaluate_models(d), "generate_pseudo_absences")
  expect_error(generate_map(d, aoi = "x.gpkg"), "generate_pseudo_absences")
})

test_that("generate_pseudo_absences demands strategy, aoi, formatted data", {
  d <- data.frame(longitude = -66.5, latitude = 47.2, year = 2023,
                  present = 1)
  expect_error(generate_pseudo_absences(d, aoi = "bbox"), "strategy")
  expect_error(generate_pseudo_absences(d, strategy = "random"), "aoi")
  expect_error(generate_pseudo_absences(d, aoi = "bbox", strategy = "sre"),
               "arg")
  raw <- data.frame(lon = -66.5, lat = 47.2, yr = 2023)
  expect_error(generate_pseudo_absences(raw, aoi = "bbox",
                                        strategy = "random"),
               "format_data")
  withabs <- rbind(d, transform(d, present = 0))
  expect_error(generate_pseudo_absences(withabs, aoi = "bbox",
                                        strategy = "random"),
               "already contains absence")
})
