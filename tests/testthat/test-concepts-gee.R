# Live Earth Engine tests. Opt in with ALPHASDM_LIVE_TESTS=1: they need
# authentication and a network, so they never run on CRAN or CI by surprise.

skip_if_no_live_gee <- function() {
  skip_if(!nzchar(Sys.getenv("ALPHASDM_LIVE_TESTS")),
          "live GEE tests are opt-in (set ALPHASDM_LIVE_TESTS=1)")
  skip_if_not(tryCatch({ ensure_gee_authenticated(); TRUE },
                       error = function(e) FALSE),
              "Earth Engine authentication unavailable")
}

test_that("the GEE encoder matches the R reference on real pixels", {
  skip_if_no_live_gee()
  ee <- reticulate::import("ee")
  sae <- make_test_sae()
  emb_cols <- sprintf("A%02d", 0:63)

  region <- ee$Geometry$Point(c(-73.6, 44.5))$buffer(5000)
  img <- get_embedding_image(2023)
  pts <- ee$FeatureCollection$randomPoints(region, 100L, 7L)
  emb_fc <- img$sampleRegions(collection = pts, scale = 10,
                              geometries = TRUE, tileScale = 16L)
  info <- retry_curl_download(emb_fc$getInfo())
  X <- do.call(rbind, lapply(info$features, function(f)
    as.numeric(unlist(f$properties[emb_cols]))))

  act_fc <- ee_concept_activations(img, sae)$
    sampleRegions(collection = emb_fc, properties = list(),
                  scale = 10, tileScale = 16L)
  act_info <- retry_curl_download(act_fc$getInfo())
  G <- do.call(rbind, lapply(act_info$features, function(f)
    as.numeric(unlist(f$properties[concept_band_names(sae$m)]))))

  R <- sae_encode(X, sae)
  # Server-side float32 against R doubles.
  expect_equal(G, unname(R), tolerance = 1e-5)
  expect_true(all(rowSums(R > 0) <= sae$k))
})

test_that("select_concepts trims the activation image to the asked-for bands", {
  skip_if_no_live_gee()
  ee <- reticulate::import("ee")
  sae <- make_test_sae()
  keep <- concept_band_names(sae$m)[c(2, 5)]
  img <- ee_concept_activations(get_embedding_image(2023), sae,
                                select_concepts = keep)
  expect_equal(unlist(img$bandNames()$getInfo()), keep)
})

test_that("derive_concepts runs end to end on a small region", {
  skip_if_no_live_gee()
  ee <- reticulate::import("ee")
  set.seed(11)
  # A handful of points near Plattsburgh NY; similarity keeps the fit light
  # and gives projection a linear direction to test.
  pts <- data.frame(longitude = stats::runif(24, -73.65, -73.55),
                    latitude  = stats::runif(24, 44.45, 44.55),
                    year = 2023, present = rep(c(1, 0), each = 12))
  aoi <- ee$Geometry$Rectangle(c(-73.7, 44.4, -73.5, 44.6))
  sae <- make_test_sae()

  fit <- suppressMessages(
    evaluate_models(pts, methods = "similarity",
                    predict_coords = pts, scale = 10))
  res <- suppressMessages(
    derive_concepts(fit, method = c("selection", "projection", "ablation"),
                    top_n = 3L, region = aoi, boot_reps = 50L,
                    sae = sae, labels = NULL, maps = TRUE))

  expect_s3_class(res, "alphasdm_concepts")
  expect_equal(nrow(res), sae$m)
  expect_true(all(c("selection_ratio", "projection_cos", "ablation_drop")
                  %in% names(res)))
  expect_equal(sum(!is.na(res$ablation_drop)), 3L)
  expect_type(attr(res, "maps"), "list")
})
