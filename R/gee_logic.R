#' Internal helper for retrying GEE operations with exponential backoff
#'
#' @noRd
retry_curl_download <- function(expr, max_retries = 5, initial_delay = 1) {
  for (i in seq_len(max_retries)) {
    res <- try(expr, silent = TRUE)
    if (!inherits(res, "try-error")) {
      return(res)
    }

    msg <- as.character(res)
    is_retryable <- grepl("429", msg) || grepl("Computation timed out", msg) ||
      grepl("Unknown Error", msg) || grepl("cannot open the connection", msg)

    if (is_retryable && i < max_retries) {
      delay <- initial_delay * (2^(i - 1)) + runif(1, 0, 1)
      sdm_warn(sprintf("GEE rate limit or timeout; retrying in %.0fs (attempt %d of %d)",
                       delay, i, max_retries), indent = 1L)
      Sys.sleep(delay)
    } else {
      stop(res)
    }
  }
}

get_embedding_image <- function(year) {
  ee <- reticulate::import("ee")
  asset_path <- "GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL"
  emb_cols <- sprintf("A%02d", 0:63)

  # No scale arguments, no reprojection. Just the raw, composited 10m image.
  # GEE will handle resampling automatically during sampling or export.
  img <- ee$ImageCollection(asset_path)$
    filter(ee$Filter$calendarRange(as.integer(year), as.integer(year), "year"))$
    mosaic()$
    select(emb_cols)

  return(img)
}

#' Sample Embeddings at FeatureCollection
#'
#' @param fc FeatureCollection with 'year' property
#' @param scale Resolution in meters
#' @param properties Optional properties to retain
#' @param geometries Boolean, retain geometries?
#' @param years Optional list of years to filter
#' @noRd
get_embeddings_at_fc_raw <- function(fc, scale, properties = NULL, geometries = FALSE, years = NULL) {
  ee <- reticulate::import("ee")

  if (is.null(years)) {
    years <- retry_curl_download(fc$aggregate_array("year")$distinct()$getInfo())
  }

  sampled_fcs <- list()
  for (yr in years) {
    yr_fc <- fc$filter(ee$Filter$eq("year", as.integer(yr)))
    yr_img <- get_embedding_image(yr)

    sampled <- yr_img$sampleRegions(
      collection = yr_fc,
      properties = as.list(properties),
      scale = scale,
      geometries = geometries,
      tileScale = 16L
    )
    sampled_fcs <- c(sampled_fcs, list(sampled))
  }
  return(ee$FeatureCollection(sampled_fcs)$flatten())
}

get_embeddings_at_fc <- function(fc, scale, properties = NULL, geometries = FALSE, years = NULL) {
  ee <- reticulate::import("ee")
  raw <- get_embeddings_at_fc_raw(fc, scale, properties, geometries, years)
  return(raw$filter(ee$Filter$notNull(as.list("A00"))))
}

#' Upload Points to GEE efficiently, chunking large DFs to stay under 10 MB
#' @param df Data frame with longitude, latitude, and optional columns.
#' @param chunk_size Max rows per GeoJSON payload (default 5000 ≈ 4 MB).
#' @return ee$FeatureCollection
upload_points_to_gee <- function(df, chunk_size = 5000L) {
  ee       <- reticulate::import("ee")
  json_mod <- reticulate::import("json")

  upload_chunk <- function(chunk_df) {
    coord_cols <- c("longitude", "latitude")
    prop_cols  <- setdiff(names(chunk_df), coord_cols)
    features   <- vector("list", nrow(chunk_df))
    for (i in seq_len(nrow(chunk_df))) {
      props <- setNames(
        lapply(prop_cols, function(col) {
          v <- chunk_df[[col]][[i]]
          if (is.integer(v)) as.integer(v) else if (is.numeric(v)) as.numeric(v) else v
        }),
        prop_cols
      )
      features[[i]] <- list(
        type       = "Feature",
        geometry   = list(type = "Point",
                          coordinates = list(chunk_df$longitude[[i]], chunk_df$latitude[[i]])),
        properties = props
      )
    }
    geojson_py <- json_mod$loads(
      as.character(jsonlite::toJSON(list(type = "FeatureCollection", features = features),
                                    auto_unbox = TRUE, digits = 10))
    )
    ee$FeatureCollection(geojson_py)
  }

  if (nrow(df) <= chunk_size) return(upload_chunk(df))

  starts <- seq(1L, nrow(df), by = chunk_size)
  fc     <- upload_chunk(df[starts[[1L]]:min(starts[[1L]] + chunk_size - 1L, nrow(df)), ])
  for (s in starts[-1L]) {
    chunk_fc <- upload_chunk(df[s:min(s + chunk_size - 1L, nrow(df)), ])
    fc       <- fc$merge(chunk_fc)
  }
  fc
}

#' Generate Background Points Server-Side
#'
#' Generates pseudo-absence points entirely on GEE. Samples only band A00
#' for land validation; it never downloads background coordinates to R.
#' Returns a GEE FeatureCollection (geometry + year + present=0).
#' @noRd
generate_background_fc_gee <- function(aoi_year, count, region) {
  ee <- reticulate::import("ee")

  sdm_info(sprintf("Target: %d background points (server-side, no pre-check)", count), indent = 1L)

  # No A00 pre-check needed: get_embeddings_at_fc already filters notNull("A00"),
  # so ocean/invalid points are dropped during the main 64-band embedding sampling.
  # This collapses two sampleRegions calls into one.
  raw_fc <- ee$FeatureCollection$randomPoints(region, as.integer(count * 3L))
  bg_fc  <- raw_fc$map(function(f) f$set("year", as.integer(aoi_year), "present", 0L))
  return(bg_fc$limit(as.integer(count)))
}

#' GEE Classifier Methods Registry
#'
#' Every supervised classifier exposed by `ee.Classifier` that can be trained
#' from labelled points is registered here. Each entry records how to build the
#' classifier and how to read its output back as a presence-suitability score
#' (higher = more suitable):
#'   fn        : the `ee.Classifier` factory name
#'   output    : value passed to `setOutputMode()`
#'   score     : band / property produced by `classify()` to read
#'   transform : how to convert the raw `score` into presence-suitability
#'     - "none"        score is already P(presence) (the SMILE probability of class 1)
#'     - "invert"      use `1 - score` (libsvm reports the probability of the
#'                     FIRST class it saw in training; see format_data row-order contract)
#'     - "mindist_raw" `score` is a RAW distance array `[d_absence, d_presence]`;
#'                     use `d_absence - d_presence` so closer-to-presence ranks higher
#'                     (minimumDistance has no PROBABILITY output mode)
#'
#' One caveat on the signed Alpha Earth embeddings: smileNaiveBayes assumes
#' positive-integer features and discards negative inputs, so it cannot use
#' embeddings that span [-1, 1]. It is registered for completeness and is not
#' a sensible choice here.
GEE_CLASSIFIER_METHODS <- list(
  rf         = list(fn = "smileRandomForest",      output = "PROBABILITY", score = "classification", transform = "none", pool = "balanced", persistable = TRUE),
  gbt        = list(fn = "smileGradientTreeBoost", output = "PROBABILITY", score = "classification", transform = "none", pool = "balanced"),
  maxent     = list(fn = "amnhMaxent",             output = "PROBABILITY", score = "probability",    transform = "none"),
  # NOTE: the svm output mode/transform below is the CLASSIFICATION-SVM fallback
  # (C_SVC / NU_SVC). When svmType is a regression SVM (the EPSILON_SVR default, or
  # NU_SVR), resolve_clf_spec() switches this to REGRESSION + transform "none" so the
  # regressed 0/1 score is read directly. See build_gee_clf_params() for the defaults.
  svm        = list(fn = "libsvm",                 output = "PROBABILITY", score = "classification", transform = "invert"),
  cart       = list(fn = "smileCart",              output = "PROBABILITY", score = "classification", transform = "none", persistable = TRUE),
  knn        = list(fn = "smileKNN",               output = "PROBABILITY", score = "classification", transform = "none", pool = "balanced"),
  naivebayes = list(fn = "smileNaiveBayes",        output = "PROBABILITY", score = "classification", transform = "none"),
  mindist    = list(fn = "minimumDistance",        output = "RAW",         score = "classification", transform = "mindist_raw")
)

#' The point pool a method trains on
#'
#' Reads the `pool` field of the registry entry, defaulting to the full pool.
#'   "balanced"  presences plus the class-balanced background. Trees would otherwise
#'               minimise loss by predicting the majority class, and smileKNN in
#'               PROBABILITY mode returns the raw positive-vote fraction, so its score
#'               IS the local class frequency: at prevalence p the expected number of
#'               positive neighbours is k * p, and once that falls below 1 the score is
#'               0 nearly everywhere.
#'   "presence"  presences only (the reducer methods, which have no negative class).
#'   "full"      presences plus all background.
#' @noRd
method_pool <- function(method) {
  if (method %in% GEE_REDUCER_METHODS) return("presence")
  pool <- GEE_CLASSIFIER_METHODS[[method]]$pool
  if (is.null(pool)) "full" else pool
}

#' Build constructor arguments for a GEE classifier
#'
#' Picks only the arguments each `ee.Classifier` factory accepts out of the
#' shared parameter list, coercing integer-typed arguments and dropping NULLs.
#' @noRd
build_gee_clf_params <- function(method, params) {
  int_or_null <- function(x) if (!is.null(x)) as.integer(x) else NULL
  # RF/GBT-oriented defaults that must not leak into other constructors.
  tree_core <- c("numberOfTrees", "minLeafPopulation", "bagFraction", "shrinkage",
                 "maxNodes", "variablesPerSplit", "lambda_", "polynomial", "batch_size")

  p <- switch(method,
    rf = list(
      numberOfTrees     = int_or_null(params$numberOfTrees),
      variablesPerSplit = int_or_null(params$variablesPerSplit),
      minLeafPopulation = int_or_null(params$minLeafPopulation),
      bagFraction       = params$bagFraction,
      maxNodes          = int_or_null(params$maxNodes)
    ),
    gbt = list(
      numberOfTrees = int_or_null(params$numberOfTrees),
      shrinkage     = params$shrinkage,
      maxNodes      = int_or_null(params$maxNodes)
    ),
    cart = list(
      maxNodes          = int_or_null(params$maxNodes),
      minLeafPopulation = int_or_null(params$minLeafPopulation)
    ),
    knn = list(
      # k sets the RESOLUTION of the output surface, not just the smoothing: PROBABILITY
      # mode returns the positive-vote fraction, so the score can only take k+1 distinct
      # values. k = 5 yields a 6-level suitability map, too coarse to rank cells or to
      # threshold sensibly; 15 keeps the neighbourhood local while giving a 16-level
      # surface. Callers with very small training sets should lower it: smile requires
      # k < n_train.
      k            = int_or_null(if (!is.null(params$k)) params$k else 15L),
      searchMethod = params$searchMethod,
      metric       = params$metric
    ),
    naivebayes = list(
      lambda = params$lambda
    ),
    mindist = list(
      metric   = params$metric,
      kNearest = int_or_null(params$kNearest)
    ),
    svm = {
      sp <- params[setdiff(names(params), tree_core)]
      # Defaults: an EPSILON_SVR with an RBF kernel (cost 10, gamma 0.05). Regressing
      # the 0/1 label gives a continuous, smoothly varying score, which suits the
      # ranking SDM needs better than the discrete class probability a C_SVC produces.
      # RBF is O(n^2); for very large training sets pass kernelType = "LINEAR" to
      # scale O(n x d).
      if (is.null(sp$svmType))    sp$svmType    <- "EPSILON_SVR"
      if (is.null(sp$kernelType)) sp$kernelType <- "RBF"
      if (is.null(sp$cost))       sp$cost       <- 10
      if (is.null(sp$gamma))      sp$gamma      <- 0.05
      # gamma is only meaningful for POLY/RBF/SIGMOID; libsvm errors if it is sent
      # with a LINEAR kernel, so drop it there.
      if (identical(sp$kernelType, "LINEAR")) sp$gamma <- NULL
      sp
    },
    maxent = params[setdiff(names(params), tree_core)],
    stop("Unsupported classifier method: ", method)
  )
  p[!vapply(p, is.null, logical(1))]
}

#' Build amnhMaxent tuning params from a regularization multiplier + feature classes
#'
#' The two ENMeval-style maxent levers (Muscarella 2014, Radosavljevic & Anderson 2014):
#' `beta` is the regularization multiplier (higher = simpler/smoother) and `features`
#' is a feature-class string: "auto" keeps GEE's sample-size-based autoFeature, else a
#' combination of L/Q/H/P/T (e.g. "LQH") turns autoFeature off and toggles those classes.
#' @noRd
maxent_tuning_params <- function(beta = 1, features = "auto") {
  mp <- list(betaMultiplier = beta)
  if (!is.null(features) && !identical(features, "auto")) {
    mp$autoFeature <- FALSE
    mp$linear      <- TRUE
    mp$quadratic   <- grepl("Q", features, ignore.case = TRUE)
    mp$hinge       <- grepl("H", features, ignore.case = TRUE)
    mp$product     <- grepl("P", features, ignore.case = TRUE)
    mp$threshold   <- grepl("T", features, ignore.case = TRUE)
  }
  mp
}

#' Resolve the output-mode / score / transform spec for a trained classifier
#'
#' Most classifiers use a static spec from `GEE_CLASSIFIER_METHODS`. `svm` is the
#' exception: libsvm can be a probabilistic classifier (C_SVC / NU_SVC) or a
#' regression (EPSILON_SVR / NU_SVR), and these are read back differently. A
#' regression SVM emits a single REGRESSION value that approximates the 0/1 label
#' (higher = more suitable), so it is read directly (transform "none") rather than
#' via the C_SVC `1 - p` probability flip.
#' @noRd
resolve_clf_spec <- function(method, filtered_params) {
  spec <- GEE_CLASSIFIER_METHODS[[method]]
  if (method == "svm") {
    svm_type <- if (!is.null(filtered_params$svmType)) filtered_params$svmType else "EPSILON_SVR"
    if (svm_type %in% c("EPSILON_SVR", "NU_SVR")) {
      spec$output    <- "REGRESSION"
      spec$score     <- "classification"
      spec$transform <- "none"
    }
  }
  spec
}

#' Methods backed by a reducer rather than an `ee.Classifier`
#' @noRd
GEE_REDUCER_METHODS <- "similarity"

#' Train GEE Model
#'
#' @param sampled_fc FeatureCollection of points carrying the A00-A63 embedding
#'   bands plus `class_property`.
#' @param method Model key: a name in `GEE_CLASSIFIER_METHODS` or `GEE_REDUCER_METHODS`.
#' @param params Named list of hyperparameters for the chosen method.
#' @param class_property Property holding the response (1 = presence, 0 = background).
#' @param persist When TRUE and the method is a persistable tree classifier
#'   (one whose registry entry sets `persistable`), the trained model is exported to a GEE
#'   asset and reloaded, the workaround for "Computed value is too large" on large
#'   random forests. The returned list then carries the `asset_id` to clean up.
#' @param project GEE project id for the temporary asset folder.
train_gee_model <- function(sampled_fc, method, params = list(), class_property = "present",
                            persist = FALSE, project = NULL) {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)

  is_classifier <- method %in% names(GEE_CLASSIFIER_METHODS)
  is_reducer <- method %in% GEE_REDUCER_METHODS

  if (!is_classifier && !is_reducer) stop("Unsupported method: ", method)

  LABEL_COL <- "label"

  sampled_fc <- sampled_fc$map(function(f) {
    f$set(LABEL_COL, ee$Number(f$get(class_property))$toInt())
  })

  # Training
  if (is_classifier) {
    clf_factory <- ee$Classifier[[GEE_CLASSIFIER_METHODS[[method]]$fn]]

    filtered_params <- build_gee_clf_params(method, params)
    spec <- resolve_clf_spec(method, filtered_params)
    clf <- do.call(clf_factory, filtered_params)
    clf <- clf$setOutputMode(spec$output)

    trained_model <- clf$train(
      features = sampled_fc,
      classProperty = LABEL_COL,
      inputProperties = emb_cols
    )

    # Optionally persist large tree models to an asset (avoids "Computed value is
    # too large" when classifying inline). Only the supported tree types qualify.
    # A reloaded classifier supports CLASSIFICATION/REGRESSION but NOT PROBABILITY, so
    # we persist in REGRESSION mode: an RF regressing the 0/1 label yields a continuous
    # suitability score (~P(presence)), which SDM ranking/AUC needs. The spec is
    # rewritten so predict reads it as a regression score.
    asset_id <- NULL
    if (persist && isTRUE(GEE_CLASSIFIER_METHODS[[method]]$persistable)) {
      # Persistence is an optimization, not a requirement: it lets a whole-region export
      # apply a STORED forest instead of retraining it on every tile. It relies on GEE's
      # batch scheduler, which under a throttled/restricted tier can be too backlogged to
      # run even this tiny export. So cap the wait short and FALL BACK to the inline
      # classifier on any failure rather than aborting the whole run; the inline
      # PROBABILITY classify path still produces a valid map (just retrains per tile).
      persisted <- tryCatch({
        # Train a FRESH regressor (you cannot re-mode a classification-trained forest),
        # export+reload it, and read it as a regression suitability score (~P(presence)).
        reg_clf <- do.call(clf_factory, filtered_params)$setOutputMode("REGRESSION")$train(
          features = sampled_fc, classProperty = LABEL_COL, inputProperties = emb_cols)
        ee_persist_classifier(reg_clf, project = project, max_minutes = 5)
      }, error = function(e) {
        sdm_warn(sprintf("Classifier persistence unavailable (%s); falling back to the inline %s classifier.",
                         conditionMessage(e), toupper(method)), indent = 1L)
        NULL
      })
      if (!is.null(persisted)) {
        trained_model <- persisted$classifier
        spec <- list(output = "REGRESSION", score = spec$score, transform = "none")
        asset_id <- persisted$asset_id
      }
    }

    return(list(
      trained       = trained_model,
      is_classifier = TRUE,
      method        = method,
      spec          = spec,
      asset_id      = asset_id
    ))
  } else {
    # Reducer logic (Simplified Centroid)
    presence_fc <- sampled_fc$filter(ee$Filter$eq(LABEL_COL, 1.0))
    res <- presence_fc$reduceColumns(
      reducer = ee$Reducer$mean()$`repeat`(64L),
      selectors = emb_cols
    )$getInfo()
    # Alpha Earth embeddings are unit-length, so a dot product between two of them is
    # the cosine of the angle between them. The mean of unit vectors is not itself
    # unit-length, though: its norm falls as the presences spread out, which would
    # scale every score for that species by an arbitrary constant. Normalising the
    # centroid makes the score a genuine cosine on [-1, 1] and comparable between
    # species. Within one species this is a constant rescale, so it leaves the
    # ranking, and therefore AUC and TSS, unchanged.
    weights <- as.numeric(res$mean)
    nrm     <- sqrt(sum(weights^2))
    if (is.finite(nrm) && nrm > 0) weights <- weights / nrm
    
    return(list(
      weights       = weights,
      is_classifier = FALSE,
      method        = method
    ))
  }
}

#' Predict GEE Map
#'
#' @param model_res Result from train_gee_model
#' @param img Alpha Earth mosaic
#' @noRd
predict_gee_map <- function(model_res, img) {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)

  if (model_res$is_classifier) {
    spec       <- if (!is.null(model_res$spec)) model_res$spec else GEE_CLASSIFIER_METHODS[[model_res$method]]
    classified <- img$classify(model_res$trained)

    if (spec$transform == "mindist_raw") {
      # RAW output is an array band [d_absence, d_presence]; presence-suitability
      # is "closer to the presence centre", i.e. d_absence - d_presence.
      flat <- classified$select(spec$score)$arrayFlatten(list(list("d_absence", "d_presence")))
      prediction <- flat$select("d_absence")$subtract(flat$select("d_presence"))
    } else {
      prediction <- classified$select(spec$score)
      if (spec$transform == "invert") {
        prediction <- ee$Image(1.0)$subtract(prediction)
      }
    }

    return(prediction$rename("similarity"))
  } else {
    # Simple Dot Product Similarity
    weights_ee <- ee$Image$constant(as.list(model_res$weights))$rename(emb_cols)
    prediction <- img$multiply(weights_ee)$reduce(ee$Reducer$sum())
    
    return(prediction$rename("similarity"))
  }
}

#' Predict Scores for Multiple Models on GEE
#'
#' @param fc GEE FeatureCollection with embeddings
#' @param models_list List of model results from AlphaSDM
#' @noRd
predict_all_models_gee <- function(fc, models_list) {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)

  methods <- names(models_list)
  classifiers <- methods[sapply(models_list, function(m) m$is_classifier)]
  reducers <- methods[!sapply(models_list, function(m) m$is_classifier)]

  scored_fc <- fc

  # Chain Classifiers
  for (m in classifiers) {
    model_res  <- models_list[[m]]
    spec       <- if (!is.null(model_res$spec)) model_res$spec else GEE_CLASSIFIER_METHODS[[model_res$method]]
    score_col  <- spec$score
    target_col <- paste0("pred_", m)
    scored_fc  <- scored_fc$classify(model_res$trained)

    if (spec$transform == "mindist_raw") {
      scored_fc <- scored_fc$map(function(f) {
        arr <- ee$Array(f$get(score_col))
        f$set(target_col, arr$get(list(0L))$subtract(arr$get(list(1L))))
      })
    } else if (spec$transform == "invert") {
      scored_fc <- scored_fc$map(function(f) f$set(target_col, ee$Number(1.0)$subtract(f$get(score_col))))
    } else {
      scored_fc <- scored_fc$map(function(f) f$set(target_col, f$get(score_col)))
    }
  }

  # Handle Reducers (Similarity)
  if (length(reducers) > 0) {
    for (m in reducers) {
        model_res <- models_list[[m]]
        centroid_ee <- ee$Array(as.list(model_res$weights))
        
        target_col <- paste0("pred_", m)
        scored_fc <- scored_fc$map(function(f) {
            point_ee <- ee$Array(f$toArray(emb_cols))
            # Linear dot product
            score <- point_ee$multiply(centroid_ee)$reduce(ee$Reducer$sum(), list(0L))$get(list(0L))
            return(f$set(target_col, score))
        })
    }
  }

  return(scored_fc)
}

#' Assign spatial cross-validation folds (blockCV blocks or k-means clusters)
#'
#' blockCV (Valavi et al. 2018) tiles the study area into square blocks (sized to a
#' fraction of the extent, or `block_size` metres) and distributes them among folds,
#' controlling spatial autocorrelation WITHOUT the harsh contiguous-region
#' extrapolation that k-means clustering of coordinates produces (k-means makes each
#' fold one geographic chunk, the most pessimistic spatial CV). Falls back to k-means
#' if `blockCV`/`sf` are unavailable or fail.
#' @noRd
assign_cv_folds <- function(df, k, method = c("block", "kmeans", "random"), block_size = NULL) {
  method <- match.arg(method)
  if (method == "random") {
    # Non-spatial k-fold: matches an interpolation (within-region) prediction task,
    # where held-out points are drawn from the same domain as training. Contrast with
    # "block", which separates folds spatially to estimate transfer/extrapolation error.
    n <- nrow(df)
    return(as.integer(sample(rep_len(seq_len(k), n), n)))
  }
  if (method == "block") {
    if (requireNamespace("blockCV", quietly = TRUE) && requireNamespace("sf", quietly = TRUE)) {
      folds <- tryCatch(blockcv_folds(df, k, block_size), error = function(e) {
        sdm_warn(sprintf("blockCV folds failed (%s); using k-means.", conditionMessage(e))); NULL })
      if (!is.null(folds) && !anyNA(folds) && length(unique(folds)) > 1) return(as.integer(folds))
    } else {
      sdm_warn("Package 'blockCV' not installed; using k-means spatial folds.")
    }
  }
  as.integer(stats::kmeans(df[, c("longitude", "latitude")], centers = k, nstart = 5L)$cluster)
}

#' blockCV spatial-block fold assignment
#' @noRd
blockcv_folds <- function(df, k, block_size = NULL) {
  sfp <- sf::st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
  ctr <- sf::st_coordinates(sf::st_centroid(sf::st_union(sfp)))
  utm <- (if (ctr[2] >= 0) 32600L else 32700L) + as.integer(floor((ctr[1] + 180) / 6)) + 1L
  sfp <- sf::st_transform(sfp, utm)                       # metric CRS for square blocks
  if (is.null(block_size)) {
    bb  <- sf::st_bbox(sfp)
    ext <- max(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
    block_size <- ext / (2 * sqrt(k) + 1)                 # several blocks per fold
  }
  sb <- blockCV::cv_spatial(x = sfp, k = as.integer(k), size = block_size,
                            selection = "random", iteration = 50L,
                            progress = FALSE, plot = FALSE, report = FALSE)
  as.integer(sb$folds_ids)
}

#' Internal: download a single-band EE image over a region as a GeoTIFF (tiled)
#'
#' Self-contained replacement for \code{rgee::ee_as_rast} that does NOT require
#' Google Drive or GCS authentication; it pulls pixels directly through the
#' synchronous \code{getDownloadURL} endpoint and mosaics the tiles locally.
#'
#' Robustness to water / missing data: AlphaEarth has no embedding over open
#' water, so predictions there are \emph{masked}. Masked pixels otherwise export
#' as ragged or zero-sized tiles, and a tile that is entirely water can come back
#' empty. We therefore \code{unmask()} the image to an explicit \code{nodata}
#' sentinel before download: every pixel (land or water) then has a real value,
#' so every tile, including all-water tiles, returns a valid, equal-sized
#' GeoTIFF. The sentinel is restored to \code{NA} in the mosaicked output.
#'
#' Tiles are sized to stay under \code{getDownloadURL}'s synchronous limits
#' (~32 MB request, 10000 px per side); a single-band float32 tile of
#' \code{max_tile_px} on a side is ~\code{max_tile_px^2 * 4} bytes.
#'
#' @param image     ee.Image with a single prediction band.
#' @param region    ee.Geometry whose bounds define the export extent.
#' @param scale     Output pixel size in metres.
#' @param dsn       Destination GeoTIFF path.
#' @param nodata    Sentinel written over masked (e.g. water) pixels, mapped to NA.
#' @param max_tile_px Tile edge length in pixels (request-size budget).
#' @param tries     Per-tile download retries (exponential backoff).
#' @return \code{dsn}; writes the GeoTIFF as a side effect.
#' @noRd
export_image_tiled <- function(image, region, scale, dsn,
                               nodata = -9999, max_tile_px = 2048L, tries = 4L) {
  ee <- reticulate::import("ee")

  # Env override for tile size: on restricted-quota projects or very large AOIs the
  # 2048 px default can exceed getDownloadURL's synchronous size/compute budget and
  # stall. ALPHASDM_MAX_TILE_PX lets a caller shrink tiles without a code change.
  .env_tile <- suppressWarnings(as.integer(Sys.getenv("ALPHASDM_MAX_TILE_PX", "")))
  if (!is.na(.env_tile) && .env_tile >= 128L) max_tile_px <- .env_tile

  # Explicit nodata so masked/water pixels download uniformly (no ragged tiles,
  # no "empty image" failures on all-water tiles). See function docs.
  img <- image$unmask(ee$Image$constant(nodata))$toFloat()

  # Region bounding box in EPSG:4326.
  ring <- region$bounds()$coordinates()$get(0L)$getInfo()
  xs   <- vapply(ring, function(p) p[[1]], numeric(1))
  ys   <- vapply(ring, function(p) p[[2]], numeric(1))
  xmin <- min(xs); xmax <- max(xs); ymin <- min(ys); ymax <- max(ys)

  # Tile size in degrees from a pixel budget (~111320 m per degree of latitude).
  tile_deg <- max_tile_px * scale / 111320
  nx <- max(1L, as.integer(ceiling((xmax - xmin) / tile_deg)))
  ny <- max(1L, as.integer(ceiling((ymax - ymin) / tile_deg)))
  xb <- seq(xmin, xmax, length.out = nx + 1L)
  yb <- seq(ymin, ymax, length.out = ny + 1L)

  # Tile scratch dir. By default a throwaway tempdir (cleaned on exit). If
  # ALPHASDM_TILE_CACHE is set, tiles stream into a PERSISTENT, per-output cache
  # keyed by the destination filename: already-downloaded, non-empty tiles are
  # skipped on a re-run, so an interrupted or throttled export resumes on just the
  # missing tiles instead of restarting. This is what makes a slow whole-region
  # export (many small tiles under a tight quota) robust to process death.
  cache_root <- Sys.getenv("ALPHASDM_TILE_CACHE", "")
  if (nzchar(cache_root)) {
    key    <- gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(basename(dsn)))
    tmpdir <- file.path(cache_root, sprintf("%s_%dm", key, as.integer(scale)))
    dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
    # persistent: do NOT unlink on exit
  } else {
    tmpdir <- file.path(tempdir(), sprintf("alphasdm_tiles_%06d", as.integer(stats::runif(1, 1, 1e6))))
    dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
    on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)
  }

  fetch_tile <- function(geom, path) {
    # resume: a cached, non-empty tile from a prior run is reused as-is
    if (file.exists(path) && file.info(path)$size > 0) return(TRUE)
    for (k in seq_len(tries)) {
      ok <- tryCatch({
        url <- img$getDownloadURL(list(region = geom, scale = scale,
                                       format = "GEO_TIFF", crs = "EPSG:4326"))
        utils::download.file(url, path, mode = "wb", quiet = TRUE)
        file.exists(path) && file.info(path)$size > 0
      }, error = function(e) FALSE)
      if (isTRUE(ok)) return(TRUE)
      Sys.sleep(2 * k)
    }
    FALSE
  }

  tiles <- character(0); failed <- 0L
  for (i in seq_len(nx)) for (j in seq_len(ny)) {
    geom <- ee$Geometry$Rectangle(c(xb[i], yb[j], xb[i + 1L], yb[j + 1L]))
    path <- file.path(tmpdir, sprintf("tile_%03d_%03d.tif", i, j))
    if (fetch_tile(geom, path)) tiles <- c(tiles, path) else failed <- failed + 1L
  }
  if (length(tiles) == 0L) stop("export_image_tiled: all tile downloads failed.")
  if (failed > 0L)
    sdm_warn(sprintf("%d of %d export tiles failed after %d retries; the output raster has gaps in those areas.",
                     failed, nx * ny, tries))

  if (length(tiles) == 1L) {
    r <- stars::read_stars(tiles[[1]], proxy = FALSE)
  } else {
    moz <- stars::st_mosaic(tiles)               # GDAL mosaic -> temp GeoTIFF
    r   <- stars::read_stars(moz, proxy = FALSE)
  }
  r[[1]][r[[1]] == nodata] <- NA                  # sentinel (water) -> NA
  stars::write_stars(r, dsn)
  dsn
}
