# The advertised workflow, verbatim (James's standing test): read a csv,
# format_data, generate_pseudo_absences, evaluate_models — only those
# calls, exactly as the README Quick Start shows them. If this test needs
# any helper code between those lines to pass, the package is broken as
# advertised, whatever the unit tests say.
# Live test: ALPHASDM_LIVE_TESTS=1 (needs Earth Engine credentials).

test_that("the README workflow runs on its own, only those lines", {
  skip_if_not(identical(Sys.getenv("ALPHASDM_LIVE_TESTS"), "1"),
              "live GEE test")

  # Stand-in for the user's file: points jittered around the Mount
  # Carleton massif (public geography, deliberately NOT survey
  # locations). The presences must actually select habitat — random
  # points ARE the regional distribution, so the environmental envelope
  # correctly contains everything and "combined" correctly refuses.
  csv <- tempfile(fileext = ".csv")
  set.seed(7)
  ridge <- rbind(c(-66.878, 47.379), c(-66.885, 47.394),
                 c(-66.900, 47.370), c(-66.860, 47.385))
  i <- sample.int(nrow(ridge), 25, replace = TRUE)
  write.csv(data.frame(lon = ridge[i, 1] + runif(25, -0.01, 0.01),
                       lat = ridge[i, 2] + runif(25, -0.008, 0.008),
                       obs_year = 2023),
            csv, row.names = FALSE)

  # ---- the advertised lines, and nothing else -------------------------------
  my_raw_df <- read.csv(csv)

  formatted_data <- format_data(
    my_raw_df,
    coords = c("lon", "lat"),
    year = "obs_year"
  )

  formatted_data <- generate_pseudo_absences(
    formatted_data,
    aoi      = list(lat = 47.38, lon = -66.88, radius = 30000),
    strategy = "combined",
    n        = 50
  )

  test_rows <- sample(nrow(formatted_data), 20)

  metrics <- evaluate_models(
    data = formatted_data[-test_rows, ],
    predict_coords = formatted_data[test_rows, ],
    scale = 10,
    aoi_year = 2023,
    methods = c("rf", "maxent")
  )
  # ---------------------------------------------------------------------------

  expect_true(all(c(1, 0) %in% formatted_data$present))
  expect_true(is.numeric(metrics$metrics$ensemble$auc_roc))
})
