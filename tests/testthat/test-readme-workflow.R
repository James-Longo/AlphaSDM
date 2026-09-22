# The advertised workflow, verbatim: the README's Example block is read out of
# README.md and run exactly as written. If it needs any helper code to pass,
# the package is broken as advertised, whatever the unit tests say.
# Live test: ALPHASDM_LIVE_TESTS=1 (needs Earth Engine credentials and network).

test_that("the README example runs as written", {
  skip_if_not(identical(Sys.getenv("ALPHASDM_LIVE_TESTS"), "1"),
              "live GEE test")
  readme <- test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README not available")

  lines <- readLines(readme)
  start <- grep("^## Example", lines)
  fence <- grep("^```", lines)
  open  <- fence[fence > start][1]
  close <- fence[fence > open][1]
  code  <- lines[(open + 1L):(close - 1L)]
  # The package under test is already loaded; library() would attach an
  # installed copy instead.
  code <- code[!grepl("^library\\(AlphaSDM\\)", code)]

  withr::local_dir(withr::local_tempdir())
  env <- new.env()
  suppressMessages(eval(parse(text = code), envir = env))

  expect_true(all(c(0, 1) %in% env$occ$present))
  expect_true(is.numeric(env$fit$metrics$ensemble$auc_roc))
  expect_true(file.exists(env$maps$ensemble_map))
})
