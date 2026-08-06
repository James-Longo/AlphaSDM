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

#' The Alpha Earth annual embedding collection
#' @noRd
ALPHAEARTH_ASSET <- "GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL"

#' Years the Alpha Earth collection covers
#'
#' Read from the collection's own image dates, so the window follows each annual
#' release instead of being pinned in the source. Cached for the session. Requires
#' a connection, like everything else in the package: there is no offline answer,
#' because a hardcoded one silently goes stale the moment Google publishes a year.
#' @noRd
alphaearth_year_range <- function() {
  if (!is.null(.alphasdm_env$year_range)) return(.alphasdm_env$year_range)
  ensure_gee_authenticated()

  ee <- reticulate::import("ee")
  ic <- ee$ImageCollection(ALPHAEARTH_ASSET)
  r  <- retry_curl_download(
    ee$List(list(ee$Date(ic$aggregate_min("system:time_start"))$get("year"),
                 ee$Date(ic$aggregate_max("system:time_start"))$get("year")))$getInfo())

  yrs <- suppressWarnings(as.integer(unlist(r)))
  if (length(yrs) != 2L || anyNA(yrs) || yrs[1] > yrs[2]) {
    stop("Could not read the Alpha Earth coverage window from ", ALPHAEARTH_ASSET, ".", call. = FALSE)
  }
  .alphasdm_env$year_range <- yrs
  yrs
}

get_embedding_image <- function(year) {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)

  # No scale argument and no reprojection: return the raw composited 10 m image.
  # Earth Engine resamples on its own when the image is sampled or exported.
  img <- ee$ImageCollection(ALPHAEARTH_ASSET)$
    filter(ee$Filter$calendarRange(as.integer(year), as.integer(year), "year"))$
    mosaic()$
    select(emb_cols)

  return(img)
}

#' Sample Embeddings at FeatureCollection
#'
#' @param fc FeatureCollection with 'year' property
#' @param scale Resolution in metres
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
#' @noRd
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
#' Generates background points entirely on Earth Engine. Samples only band A00
#' for land validation; it never downloads background coordinates to R.
#' Returns a GEE FeatureCollection (geometry + year + present=0).
#' @noRd
generate_background_fc_gee <- function(aoi_year, count, region) {
  ee <- reticulate::import("ee")

  sdm_info(sprintf("Target: %d background points (server-side, no pre-check)", count), indent = 1L)

  # No A00 pre-check here. get_embeddings_at_fc() already filters notNull("A00"),
  # so points over water and other gaps drop out during the main 64-band sampling.
  # Checking first would mean two sampleRegions calls where one does.
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
#' @noRd
GEE_CLASSIFIER_METHODS <- list(
  rf         = list(fn = "smileRandomForest",      output = "PROBABILITY", score = "classification", transform = "none", pool = "balanced", persistable = TRUE),
  gbt        = list(fn = "smileGradientTreeBoost", output = "PROBABILITY", score = "classification", transform = "none", pool = "balanced"),
  maxent     = list(fn = "amnhMaxent",             output = "PROBABILITY", score = "probability",    transform = "none"),
  # The svm entry below is the classification-SVM case, C_SVC or NU_SVC. For a
  # regression SVM, which is the EPSILON_SVR default and also NU_SVR,
  # resolve_clf_spec() swaps it to REGRESSION with transform "none" so the regressed
  # 0/1 score is read directly. build_gee_clf_params() holds the defaults.
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
  # Tree arguments. libsvm and amnhMaxent reject these, so they are stripped for
  # those two rather than listed per method.
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
      # k sets the resolution of the output surface, not just its smoothness. In
      # PROBABILITY mode the score is the positive-vote fraction among k neighbours,
      # so it can take only k + 1 distinct values. At k = 5 that is a 6-level
      # suitability map, too coarse to rank cells or to threshold. 15 keeps the
      # neighbourhood local and gives 16 levels. Lower it for a small training set:
      # smile requires k < n_train.
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
      # Defaults: EPSILON_SVR with an RBF kernel, cost 10, gamma 0.05. Regressing the
      # 0/1 label gives a continuous score, which suits ranking better than the
      # discrete class probability a C_SVC produces. RBF is O(n^2) in the number of
      # training points, so for a very large training set pass kernelType = "LINEAR",
      # which scales as O(n x d).
      if (is.null(sp$svmType))    sp$svmType    <- "EPSILON_SVR"
      if (is.null(sp$kernelType)) sp$kernelType <- "RBF"
      if (is.null(sp$cost))       sp$cost       <- 10
      if (is.null(sp$gamma))      sp$gamma      <- 0.05
      # gamma applies only to the POLY, RBF and SIGMOID kernels. libsvm errors if it
      # is sent with a LINEAR kernel, so drop it there.
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
#' @noRd
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

    # Optionally store a large tree model in an asset, which avoids "Computed value
    # is too large" when classifying inline. Only the tree types marked persistable
    # in the registry qualify.
    #
    # A reloaded classifier supports CLASSIFICATION and REGRESSION but not
    # PROBABILITY, so store it in REGRESSION mode. A forest regressing the 0/1 label
    # returns a continuous score close to the probability of presence, which is what
    # ranking needs. The spec is rewritten so predict reads it as a regression score.
    asset_id <- NULL
    if (persist && isTRUE(GEE_CLASSIFIER_METHODS[[method]]$persistable)) {
      # Storing the model is an optimisation, not a requirement. It lets a
      # whole-region export apply a stored forest instead of retraining it for every
      # tile, and it keeps the classify graph small enough that the interactive
      # endpoint will serve a tile at all. It depends on the Earth Engine batch
      # scheduler, so fall back to the inline classifier on any failure: that path
      # still produces a valid map.
      persisted <- tryCatch({
        # Train a fresh regressor: a forest already trained for classification cannot
        # be switched to another output mode. Export it, load it back, and read it as
        # a regression score.
        reg_clf <- do.call(clf_factory, filtered_params)$setOutputMode("REGRESSION")$train(
          features = sampled_fc, classProperty = LABEL_COL, inputProperties = emb_cols)
        # No limit. Every number tried here was a guess, and each one eventually cut
        # off a task that was working: a large training set takes longer to write than
        # a small one, and the write is what keeps the classify graph small enough to
        # serve a tile at all. Losing it is what forces the expensive path.
        ee_persist_classifier(reg_clf, project = project)
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
      # RAW output is an array band holding [d_absence, d_presence]. Suitability
      # means closer to the presence centre, so take d_absence - d_presence.
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

  if (length(reducers) > 0) {
    for (m in reducers) {
        model_res <- models_list[[m]]
        centroid_ee <- ee$Array(as.list(model_res$weights))
        
        target_col <- paste0("pred_", m)
        scored_fc <- scored_fc$map(function(f) {
            point_ee <- ee$Array(f$toArray(emb_cols))
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
    # Non-spatial k-fold suits an interpolation task, where the points to be
    # predicted come from the same region as the training points. Use "block"
    # instead to estimate error when transferring to a region not sampled.
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
  MIN_TILE_PX <- 128L

  # No download deadline. R's default is 60 seconds, and options(timeout) is
  # libcurl's total transfer time rather than an idle timeout, so any value at all
  # eventually cuts off a download that is healthy but slow. A tile at fine scale
  # takes as long as Earth Engine takes. Let the server decide when a request is
  # finished or refused; 0 removes the limit.
  old_timeout <- getOption("timeout")
  options(timeout = 0)
  on.exit(options(timeout = old_timeout), add = TRUE)

  # Set ALPHASDM_MAX_TILE_PX to pin the tile size. Otherwise the size is found by
  # trying and halving, below.
  pinned <- suppressWarnings(as.integer(Sys.getenv("ALPHASDM_MAX_TILE_PX", "")))
  pinned <- !is.na(pinned) && pinned >= MIN_TILE_PX
  if (pinned) max_tile_px <- suppressWarnings(as.integer(Sys.getenv("ALPHASDM_MAX_TILE_PX")))

  # Write an explicit nodata value so masked pixels download like any other. Without
  # it tiles come back ragged, and a tile that is entirely water fails as an empty
  # image. See the function documentation above.
  img <- image$unmask(ee$Image$constant(nodata))$toFloat()

  # Bounding box of the region, in EPSG:4326.
  ring <- region$bounds()$coordinates()$get(0L)$getInfo()
  xs   <- vapply(ring, function(p) p[[1]], numeric(1))
  ys   <- vapply(ring, function(p) p[[2]], numeric(1))
  xmin <- min(xs); xmax <- max(xs); ymin <- min(ys); ymax <- max(ys)

  # Tile grid for a given tile size, at about 111320 m per degree of latitude.
  build_grid <- function(px) {
    tile_deg <- px * scale / 111320
    nx <- max(1L, as.integer(ceiling((xmax - xmin) / tile_deg)))
    ny <- max(1L, as.integer(ceiling((ymax - ymin) / tile_deg)))
    list(nx = nx, ny = ny,
         xb = seq(xmin, xmax, length.out = nx + 1L),
         yb = seq(ymin, ymax, length.out = ny + 1L))
  }

  # Scratch directory for tiles. The default is a temporary directory, removed on
  # exit. Set ALPHASDM_TILE_CACHE to keep the tiles instead, in a cache keyed by the
  # destination, the scale and the tile size. A re-run then skips tiles already
  # downloaded, so an export interrupted partway resumes on the missing tiles alone.
  # That is what lets a slow whole-region export survive the process being killed.
  # The tile size is part of the key because changing it changes the grid, which
  # would otherwise leave tiles from a different grid in place.
  cache_root <- Sys.getenv("ALPHASDM_TILE_CACHE", "")
  tile_dir <- function(px) {
    if (nzchar(cache_root)) {
      key <- gsub("[^A-Za-z0-9]+", "_", tools::file_path_sans_ext(basename(dsn)))
      d <- file.path(cache_root, sprintf("%s_%dm_%dpx", key, as.integer(scale), px))
      dir.create(d, showWarnings = FALSE, recursive = TRUE)
      d
    } else {
      d <- file.path(tempdir(), sprintf("alphasdm_tiles_%06d_%dpx",
                                        as.integer(stats::runif(1, 1, 1e6)), px))
      dir.create(d, showWarnings = FALSE, recursive = TRUE)
      d
    }
  }

  # Fetch one tile. Returns TRUE, or the condition message so the caller can tell a
  # server rejection from a transient network failure.
  # A tile that cannot be served is one the batch system should compute instead.
  # is_gee_timeout() already lists the messages Earth Engine uses for a request that
  # was too expensive, including "Computation timed out" and "User memory limit
  # exceeded". A 400 is the same refusal arriving as an HTTP status, and "Timeout of
  # N seconds" is R's own download giving up. None is fixed by asking again.
  hopeless <- function(msg) {
    is_gee_timeout(msg) ||
      grepl("400", msg, fixed = TRUE) ||
      grepl("Timeout of", msg, fixed = TRUE)
  }

  fetch_tile <- function(geom, path, attempts, stop_when_hopeless = FALSE) {
    if (file.exists(path) && file.info(path)$size > 0) return(TRUE)
    last <- "unknown error"
    for (k in seq_len(attempts)) {
      # download.file reports an HTTP status as a warning and only sometimes as an
      # error, so catch the warning text too. Discarding it leaves nothing to tell a
      # server refusal from an ordinary empty result, which is the difference between
      # shrinking the tile and giving up.
      warn <- ""
      res <- tryCatch(
        withCallingHandlers({
          url <- img$getDownloadURL(list(region = geom, scale = scale,
                                         format = "GEO_TIFF", crs = "EPSG:4326"))
          utils::download.file(url, path, mode = "wb", quiet = TRUE)
          if (file.exists(path) && file.info(path)$size > 0) TRUE else "empty response"
        }, warning = function(w) {
          warn <<- conditionMessage(w); invokeRestart("muffleWarning")
        }),
        error = function(e) conditionMessage(e))
      if (isTRUE(res)) return(TRUE)
      if (nzchar(warn)) res <- warn
      last <- res
      if (hopeless(last) && stop_when_hopeless) return(last)
      if (grepl("400", last, fixed = TRUE)) return(last)
      Sys.sleep(2 * k)
    }
    last
  }

  # Try the first tile. A refusal means the interactive endpoint will not serve this
  # image at this tile size, which at 10 m with a fitted model is the normal case,
  # since every tile download re-runs the model over that tile's pixels.
  #
  # When that happens, compute the image once through Earth Engine's batch system,
  # which has a larger memory allowance and a longer timeout, and read tiles from the
  # stored result instead. Reading stored pixels costs a read, so the tile size that
  # failed will now succeed. This is Earth Engine's own advice for "User memory limit
  # exceeded": run it as an Export and read back the exported image.
  #
  # If even the stored image will not serve at this tile size, fall back to halving
  # the tile, which is the older behaviour and still correct, just slower.
  grid <- NULL; tmpdir <- NULL; first_path <- NULL; switched <- FALSE
  # The staged copy is deleted as soon as every tile is on disk. This keeps it as a
  # backstop for an early exit only.
  staged <- NULL
  on.exit(if (!is.null(staged)) ee_delete_asset_quietly(staged), add = TRUE)
  repeat {
    grid   <- build_grid(max_tile_px)
    tmpdir <- tile_dir(max_tile_px)
    first_path <- file.path(tmpdir, "tile_001_001.tif")
    geom <- ee$Geometry$Rectangle(c(grid$xb[1], grid$yb[1], grid$xb[2], grid$yb[2]))
    # fetch_tile returns a 400 immediately without spending retries, so a full retry
    # budget here still separates a server refusal from a transient network failure.
    res  <- fetch_tile(geom, first_path, attempts = tries, stop_when_hopeless = TRUE)
    if (isTRUE(res)) break
    # Shrink only when the server refused the request. A slow tile is not a tile that
    # is too large, and treating it as one used to quadruple the tile count over a
    # whole region on the strength of one slow download.
    too_big <- hopeless(res)
    if (!too_big) {
      stop(sprintf("export_image_tiled: tile download failed at %d px. %s",
                   max_tile_px, res), call. = FALSE)
    }

    # First refusal: move the computation into the batch system and try again from
    # the stored image, at the same tile size.
    if (!switched) {
      switched <- TRUE
      sdm_warn(paste("This map is too heavy to download directly at",
                     "this resolution. Computing it through Earth Engine's batch",
                     "system instead, then downloading the result."), indent = 2L)
      sdm_info("This runs server-side and can take a while. Progress is reported below.",
               indent = 2L)
      stored <- tryCatch(ee_persist_image(img, region, scale, project = NULL),
                         error = function(e) e)
      if (inherits(stored, "error")) {
        sdm_warn(sprintf("Batch export unavailable (%s); falling back to smaller tiles.",
                         conditionMessage(stored)), indent = 2L)
      } else {
        img <- stored$img
        staged <- stored$asset_id
        sdm_info("Map computed. Downloading it now.", indent = 2L)
        next                       # retry the probe against the stored image
      }
    }

    if (pinned || max_tile_px <= MIN_TILE_PX) {
      stop(sprintf("export_image_tiled: tile download failed at %d px. %s",
                   max_tile_px, res), call. = FALSE)
    }
    max_tile_px <- max(MIN_TILE_PX, max_tile_px %/% 2L)
    sdm_info(sprintf("Tile of %d px did not render; retrying at %d px.",
                     max_tile_px * 2L, max_tile_px), indent = 2L)
  }
  if (!nzchar(cache_root)) on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  nx <- grid$nx; ny <- grid$ny; xb <- grid$xb; yb <- grid$yb
  n_total <- nx * ny
  if (n_total > 50L) {
    sdm_info(sprintf("Exporting %d tiles at %d px (%.1f km each).",
                     n_total, max_tile_px, max_tile_px * scale / 1000), indent = 2L)
  }

  # Fetch every tile. Already-downloaded tiles are reused, so this is safe to call a
  # second time after switching to the stored image: only the missing tiles are
  # fetched again.
  run_tiles <- function() {
    tiles <- character(0); failed <- 0L; last_msg <- ""; refused <- FALSE
    t0 <- proc.time()[["elapsed"]]; done <- 0L
    for (i in seq_len(nx)) for (j in seq_len(ny)) {
      geom <- ee$Geometry$Rectangle(c(xb[i], yb[j], xb[i + 1L], yb[j + 1L]))
      path <- file.path(tmpdir, sprintf("tile_%03d_%03d.tif", i, j))
      res  <- fetch_tile(geom, path, attempts = tries)
      if (isTRUE(res)) {
        tiles <- c(tiles, path)
      } else {
        failed <- failed + 1L; last_msg <- res
        if (hopeless(res)) refused <- TRUE
      }
      done <- done + 1L
      # Report progress on a long export, with an estimate from the rate so far.
      if (n_total > 50L && done %% max(10L, n_total %/% 20L) == 0L) {
        el <- proc.time()[["elapsed"]] - t0
        sdm_info(sprintf("%d/%d tiles (%.0f%%), %.0f min elapsed, about %.0f min left",
                         done, n_total, 100 * done / n_total, el / 60,
                         (el / done) * (n_total - done) / 60), indent = 2L)
      }
    }
    list(tiles = tiles, failed = failed, last_msg = last_msg, refused = refused)
  }

  r <- run_tiles()

  # A tile can be refused partway through even when the first one was served: a
  # busier region, or the tier throttling as the export goes on. Escalate here too,
  # rather than leaving holes in the raster. The tiles already on disk are kept, so
  # only the missing ones are fetched from the stored image.
  if (r$refused && !switched) {
    switched <- TRUE
    sdm_warn(paste("Tiles started being refused partway through. Computing the map",
                   "through Earth Engine's batch system instead, then finishing the",
                   "download."), indent = 2L)
    stored <- tryCatch(ee_persist_image(img, region, scale, project = NULL),
                       error = function(e) e)
    if (inherits(stored, "error")) {
      sdm_warn(sprintf("Batch export unavailable (%s); keeping what downloaded.",
                       conditionMessage(stored)), indent = 2L)
    } else {
      img <- stored$img
      staged <- stored$asset_id
      sdm_info(sprintf("Map computed. Fetching the %d remaining tile%s.",
                       r$failed, if (r$failed == 1L) "" else "s"), indent = 2L)
      r <- run_tiles()
    }
  }
  tiles <- r$tiles; failed <- r$failed; last_msg <- r$last_msg

  # Every tile that is going to arrive has arrived, so the staged copy has done its
  # job. Release it now rather than at function exit: mosaicking and writing a large
  # raster takes time, and there is no reason to hold Earth Engine storage during it.
  if (!is.null(staged) && length(tiles) > 0L) {
    on_disk <- sum(file.exists(tiles) & file.info(tiles)$size > 0)
    if (on_disk == length(tiles)) {
      ee_delete_asset_quietly(staged)
      staged <- NULL
      sdm_info(sprintf("%d tile%s confirmed on disk; removed the temporary copy from your Earth Engine storage.",
                       on_disk, if (on_disk == 1L) "" else "s"), indent = 2L)
    }
  }

  if (length(tiles) == 0L)
    stop(sprintf("export_image_tiled: every tile failed. Last error: %s", last_msg), call. = FALSE)
  if (failed > 0L)
    sdm_warn(sprintf("%d of %d export tiles failed after %d retries; the raster has gaps there. Last error: %s",
                     failed, n_total, tries, last_msg))

  # Deliver the tiles rather than one mosaicked raster. A whole-region raster at fine
  # scale does not fit in memory: the Greater Antilles at 10 m is 62 GB, and reading
  # that to substitute a nodata value would fail on any ordinary machine. Each tile is
  # bounded, so the same work per tile always fits.
  #
  # The sentinel written for masked pixels is turned into NA here, one tile at a time,
  # so the tiles are usable as they are. Mosaic them with gdalbuildvrt or terra::vrt
  # if a single raster is wanted.
  out_dir <- file.path(dirname(dsn), paste0(tools::file_path_sans_ext(basename(dsn)), "_tiles"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  written <- character(0)
  for (k in seq_along(tiles)) {
    dest <- file.path(out_dir, basename(tiles[k]))
    ok <- tryCatch({
      r <- stars::read_stars(tiles[k], proxy = FALSE)
      r[[1]][r[[1]] == nodata] <- NA
      stars::write_stars(r, dest)
      TRUE
    }, error = function(e) FALSE)
    if (isTRUE(ok)) {
      written <- c(written, dest)
      # Release the downloaded copy as soon as it has been converted. Holding both
      # doubles peak disk, which at fine scale over a large region is tens of GB.
      # Tiles in the resumable cache are kept, since that is what makes a killed run
      # restartable.
      if (!nzchar(cache_root)) unlink(tiles[k])
    }
    if (n_total > 50L && k %% max(10L, length(tiles) %/% 10L) == 0L)
      sdm_info(sprintf("wrote %d/%d tiles", k, length(tiles)), indent = 2L)
  }
  if (length(written) == 0L) stop("export_image_tiled: no tile could be written.", call. = FALSE)

  # A single tile is the whole region, so hand back the raster itself rather than a
  # directory holding one file.
  if (length(written) == 1L && n_total == 1L) {
    file.copy(written[1], dsn, overwrite = TRUE)
    unlink(out_dir, recursive = TRUE)
    return(dsn)
  }
  sdm_info(sprintf("%d tiles written to %s", length(written), out_dir), indent = 2L)
  out_dir
}
