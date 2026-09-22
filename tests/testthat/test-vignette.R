# The AlphaSDM vignette, run as written. It is precomputed (the
# knitted AlphaSDM.Rmd is committed because CRAN cannot reach Earth
# Engine), so this runs the source, AlphaSDM.Rmd.orig, to keep the
# documented example from going stale.
# Live test: ALPHASDM_LIVE_TESTS=1 (needs Earth Engine credentials and network).

test_that("the AlphaSDM vignette runs end to end", {
  skip_if_not(identical(Sys.getenv("ALPHASDM_LIVE_TESTS"), "1"),
              "live GEE test")
  skip_if_not_installed("knitr")
  rmd <- test_path("..", "..", "vignettes", "AlphaSDM.Rmd.orig")
  skip_if_not(file.exists(rmd), "vignette source not available")

  script <- withr::local_tempfile(fileext = ".R")
  knitr::purl(rmd, output = script, quiet = TRUE, documentation = 0)
  code <- readLines(script)
  # The package under test is already loaded; library() would attach an
  # installed copy instead.
  code <- code[!grepl("^library\\(AlphaSDM\\)", code)]

  withr::local_pdf(withr::local_tempfile(fileext = ".pdf"))
  env <- new.env()
  suppressMessages(eval(parse(text = code), envir = env))

  expect_gt(sum(env$occ$present == 1), 0)
  expect_gt(sum(env$occ$present == 0), 0)
  expect_true("ensemble" %in% names(env$fit$metrics))
  expect_true(file.exists(env$maps$ensemble_map))
})
