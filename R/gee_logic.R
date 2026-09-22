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

#' The 64 embedding band names, A00 to A63
#' @noRd
EMB_BANDS <- sprintf("A%02d", 0:63)

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
  # No scale argument and no reprojection: return the raw composited 10 m image.
  # Earth Engine resamples on its own when the image is sampled or exported.
  img <- ee$ImageCollection(ALPHAEARTH_ASSET)$
    filter(ee$Filter$calendarRange(as.integer(year), as.integer(year), "year"))$
    mosaic()$
    select(EMB_BANDS)

  return(img)
}

#' Sample the embeddings at a FeatureCollection
#'
#' Points whose pixel is masked (open water, no coverage) are dropped.
#' @param fc FeatureCollection with a 'year' property
#' @param scale Resolution in metres
#' @param properties Optional properties to retain
#' @param geometries Boolean, retain geometries?
#' @param years List of the years present in `fc`
#' @noRd
get_embeddings_at_fc <- function(fc, scale, properties = NULL, geometries = FALSE, years) {
  ee <- reticulate::import("ee")
  sampled_fcs <- lapply(years, function(yr) {
    get_embedding_image(yr)$sampleRegions(
      collection = fc$filter(ee$Filter$eq("year", as.integer(yr))),
      properties = as.list(properties),
      scale = scale,
      geometries = geometries,
      tileScale = 16L
    )
  })
  ee$FeatureCollection(sampled_fcs)$flatten()$filter(ee$Filter$notNull(list("A00")))
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
    ee$FeatureCollection(json_mod$loads(points_geojson(chunk_df)))
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

#' GeoJSON FeatureCollection text for a point data frame
#'
#' Built column-wise: one Point per row from `longitude`/`latitude`, every
#' other column a property. Integers stay integers; doubles keep 15
#' significant digits.
#' @noRd
points_geojson <- function(df) {
  json_vals <- function(v) {
    out <- if (is.integer(v)) as.character(v)
           else if (is.numeric(v)) sprintf("%.15g", v)
           else vapply(as.character(v), function(x)
             as.character(jsonlite::toJSON(x, auto_unbox = TRUE)), character(1))
    out[is.na(v)] <- "null"
    out
  }
  prop_cols <- setdiff(names(df), c("longitude", "latitude"))
  props <- if (length(prop_cols))
    do.call(paste, c(lapply(prop_cols, function(col)
      paste0("\"", col, "\":", json_vals(df[[col]]))), sep = ","))
  else rep("", nrow(df))
  feats <- sprintf(
    '{"type":"Feature","geometry":{"type":"Point","coordinates":[%s,%s]},"properties":{%s}}',
    json_vals(df$longitude), json_vals(df$latitude), props)
  paste0('{"type":"FeatureCollection","features":[', paste(feats, collapse = ","), "]}")
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
  # glm follows the paper's regression recipe (large random pool, equal
  # total class weights), so it trains on the full pool like maxent.
  if (identical(method, "glm")) return("full")
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
GEE_REDUCER_METHODS <- c("similarity", "glm")

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
  emb_cols <- EMB_BANDS

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
  } else if (identical(method, "glm")) {
    # Logistic GLM fitted by IRLS, entirely server-side: each iteration is
    # one weighted least-squares solve (ee.Reducer.linearRegression) over
    # the training table, with the working response and weights recomputed
    # from the current coefficients. Fixed iteration count keeps the
    # round-trips deterministic. Class weights equalize the total presence
    # and absence weight (Barbet-Massin et al. 2012, GLM recipe).
    counts <- sampled_fc$aggregate_histogram(LABEL_COL)$getInfo()
    n1 <- as.numeric(counts[["1"]]); n0 <- as.numeric(counts[["0"]])
    # The only enforced bound is the mathematical one: 65 coefficients
    # (intercept + 64 bands) need more rows than coefficients or the
    # least-squares solve is undefined. How many MORE is the user's
    # judgment (the literature recipe is ~10,000 random pseudo-absences).
    if (is.na(n1) || is.na(n0) || n1 + n0 <= 65)
      stop("glm fits 65 coefficients (intercept + 64 embedding bands) and ",
           "has only ", n1 + n0, " training rows; the solve is undefined ",
           "with rows <= coefficients.", call. = FALSE)
    cw0 <- n1 / n0
    xcols <- sprintf("X%02d", 0:64)
    beta <- rep(0, 65)
    for (it in seq_len(8L)) {
      b0 <- beta[1]
      bw <- ee$Array(as.list(unname(beta[-1])))
      prepared <- sampled_fc$map(function(f) {
        xa  <- ee$Array(f$toArray(emb_cols))
        eta <- ee$Number(xa$multiply(bw)$reduce(ee$Reducer$sum(), list(0L))$
                           get(list(0L)))$add(b0)
        mu  <- ee$Number(1)$divide(eta$multiply(-1)$exp()$add(1))$
          clamp(1e-6, 1 - 1e-6)
        y   <- ee$Number(f$get(LABEL_COL))
        cw  <- y$multiply(1 - cw0)$add(cw0)      # 1 for presence, cw0 for bg
        wgt <- mu$multiply(ee$Number(1)$subtract(mu))$multiply(cw)$max(1e-6)
        z   <- eta$add(y$subtract(mu)$divide(wgt$divide(cw)))
        sw  <- wgt$sqrt()
        xs  <- ee$List(list(sw))$cat(xa$multiply(sw)$toList())
        f$set(ee$Dictionary$fromLists(as.list(xcols), xs))$
          set("Yw", z$multiply(sw))
      })
      fit <- prepared$reduceColumns(
        reducer = ee$Reducer$linearRegression(65L, 1L),
        selectors = as.list(c(xcols, "Yw")))$get("coefficients")
      beta <- as.numeric(unlist(retry_curl_download(fit$getInfo())))
      if (any(!is.finite(beta)))
        stop("glm IRLS diverged (non-finite coefficients); the classes may ",
             "be perfectly separable in embedding space.", call. = FALSE)
    }
    return(list(
      weights       = beta[-1],
      intercept     = beta[1],
      link          = "logistic",
      is_classifier = FALSE,
      method        = method
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
  emb_cols <- EMB_BANDS

  if (model_res$is_classifier) {
    spec       <- model_res$spec
    if (!is.null(model_res$replicates) && spec$transform == "none") {
      # Replicate models: classify once per model, average the maps.
      imgs <- lapply(c(list(model_res$trained), model_res$replicates),
                     function(tr) img$classify(tr)$select(spec$score))
      return(ee$ImageCollection(imgs)$mean()$rename("similarity"))
    }
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
    if (!is.null(model_res$intercept))
      prediction <- prediction$add(model_res$intercept)
    if (identical(model_res$link, "logistic"))
      prediction <- ee$Image(1)$divide(
        prediction$multiply(-1)$exp()$add(1))

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
  emb_cols <- EMB_BANDS

  methods <- names(models_list)
  classifiers <- methods[sapply(models_list, function(m) m$is_classifier)]
  reducers <- methods[!sapply(models_list, function(m) m$is_classifier)]

  scored_fc <- fc

  for (m in classifiers) {
    model_res  <- models_list[[m]]
    spec       <- model_res$spec
    score_col  <- spec$score
    target_col <- paste0("pred_", m)

    # Replicate models (balanced-pool methods trained on k background
    # draws) score each feature once per model; the prediction is the
    # mean. Only transform "none" methods carry replicates. Each score is
    # copied into its own property immediately: classify() must not be
    # given an outputName, because a LOADED classifier ignores it and
    # writes "classification" regardless, and successive classifies
    # overwrite it. A persisted main model scores in REGRESSION mode and
    # its replicates in PROBABILITY; both approximate P(presence).
    if (!is.null(model_res$replicates) && spec$transform == "none") {
      all_tr <- c(list(model_res$trained), model_res$replicates)
      rprops <- paste0("rep_", seq_along(all_tr))
      reg_score <- GEE_CLASSIFIER_METHODS[[model_res$method]]$score
      for (i in seq_along(all_tr)) {
        prop_i <- if (i == 1L) score_col else reg_score
        rp_i <- rprops[i]
        scored_fc <- scored_fc$classify(all_tr[[i]])$map(
          local({
            p <- prop_i; rp <- rp_i
            function(f) f$set(rp, f$get(p))
          }))
      }
      k_ee <- length(all_tr)
      scored_fc <- scored_fc$map(function(f) {
        tot <- ee$Number(0)
        for (rp in rprops) tot <- tot$add(ee$Number(f$get(rp)))
        f$set(target_col, tot$divide(k_ee))
      })
      next
    }
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
        b0 <- model_res$intercept
        logistic <- identical(model_res$link, "logistic")

        target_col <- paste0("pred_", m)
        scored_fc <- scored_fc$map(function(f) {
            point_ee <- ee$Array(f$toArray(emb_cols))
            score <- ee$Number(point_ee$multiply(centroid_ee)$
              reduce(ee$Reducer$sum(), list(0L))$get(list(0L)))
            if (!is.null(b0)) score <- score$add(b0)
            if (logistic)
              score <- ee$Number(1)$divide(score$multiply(-1)$exp()$add(1))
            return(f$set(target_col, score))
        })
    }
  }

  return(scored_fc)
}


#' Messages that mean a map tile was too expensive to serve
#'
#' A refusal is not fixed by asking again, only by asking for less (a smaller
#' tile) or by moving the work to the batch system. Besides Earth Engine's own
#' compute-limit messages, a 400 is the same refusal arriving as an HTTP status,
#' and "Total request size" is a tile too large to download at all.
#' @noRd
GEE_REFUSAL_PATTERN <- paste(GEE_LIMIT_PATTERN, "400", "Total request size",
                             "Timeout of", sep = "|")

#' Whether an Earth Engine error means the request was too expensive
#' @noRd
ee_refused <- function(msg) grepl(GEE_REFUSAL_PATTERN, msg, ignore.case = TRUE)

#' Download one map tile, retrying transient failures
#'
#' Runs in a worker process of export_image()'s download pool, so it is
#' self-contained: base R only, and everything it needs arrives in `job`
#' (`url`, `path`, `tries`, `refusal_pattern`).
#' @return TRUE, or the last error message.
#' @noRd
fetch_tile <- function(job) {
  options(timeout = 0)  # a tile takes as long as Earth Engine takes to compute it
  if (!is.null(job$error)) return(job$error)
  last <- "unknown error"
  for (k in seq_len(job$tries)) {
    warn <- ""
    res <- tryCatch(withCallingHandlers({
      utils::download.file(job$url, job$path, mode = "wb", quiet = TRUE)
      ok <- file.exists(job$path) && file.size(job$path) > 8 &&
        rawToChar(readBin(job$path, "raw", 2L)) %in% c("II", "MM")
      if (ok) TRUE else "empty response"
    }, warning = function(w) {
      warn <<- conditionMessage(w); invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
    if (isTRUE(res)) return(TRUE)
    last <- if (nzchar(warn)) warn else res
    if (grepl(job$refusal_pattern, last, ignore.case = TRUE)) return(last)
    Sys.sleep(2^k)
  }
  last
}

#' Download a multi-band image as one GeoTIFF per band
#'
#' The direct route. The region is cut into tiles on one pixel grid, each as
#' large as a single `getDownloadURL` request allows (32 MB), requested
#' several at a time, and stitched locally with GDAL. A tile Earth Engine
#' refuses as too expensive (memory or compute, which a fitted model can hit
#' well below the size limit) is split into quarters and retried, and tiles
#' still waiting are shrunk to match. Only a tile refused even at the minimum
#' size sends the region to Earth Engine's batch system, which writes tiles to
#' Google Drive (`ee_export_image_drive()`); Earth Engine recommends Export for
#' work too large for interactive requests.
#'
#' AlphaEarth is masked over open water, so masked pixels are written as the
#' `nodata` value and flagged as nodata in the output files.
#'
#' @param image  ee.Image; its bands, in order, are written to `dsn`.
#' @param region ee.Geometry whose bounds define the export extent.
#' @param scale  Output pixel size in metres.
#' @param dsn    Character vector of output GeoTIFF paths, one per band.
#' @param nodata Value written for masked pixels.
#' @return `dsn`.
#' @noRd
export_image <- function(image, region, scale, dsn, nodata = -9999) {
  ee <- reticulate::import("ee")
  MAX_BYTES   <- 30e6  # under getDownloadURL's 32 MB response cap
  MIN_TILE_PX <- 128L
  # Requests in flight. Measured on a full 10 m map: 16 was no faster than 8,
  # since Earth Engine limits each account's compute, and 8 stays well under the
  # concurrent-request quota.
  CONCURRENT  <- 8L
  TRIES       <- 5L

  # fetch_tile() lifts the download deadline; restore the caller's afterwards.
  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)

  # One global pixel grid in EPSG:4326, so tiles line up exactly when stitched.
  ring <- region$bounds()$coordinates()$get(0L)$getInfo()
  xs  <- vapply(ring, function(p) p[[1]], numeric(1))
  ys  <- vapply(ring, function(p) p[[2]], numeric(1))
  dpp <- scale / 111320
  x0  <- min(xs); y0 <- max(ys)
  W   <- max(1L, as.integer(ceiling((max(xs) - x0) / dpp)))
  H   <- max(1L, as.integer(ceiling((y0 - min(ys)) / dpp)))

  n_bands  <- length(dsn)
  cap_px   <- max(MIN_TILE_PX, as.integer(floor(sqrt(MAX_BYTES / (4 * n_bands)))))
  img      <- image$unmask(nodata, FALSE)$toFloat()
  tile_dir <- tempfile("alphasdm_tiles_")
  dir.create(tile_dir)
  on.exit(unlink(tile_dir, recursive = TRUE), add = TRUE)

  # A tile is its pixel offset and size on the grid.
  cut <- function(c0, r0, w, h, px) {
    cs <- seq(c0, c0 + w - 1L, by = px); rs <- seq(r0, r0 + h - 1L, by = px)
    unlist(lapply(cs, function(c) lapply(rs, function(r)
      c(c, r, min(px, c0 + w - c), min(px, r0 + h - r)))), recursive = FALSE)
  }
  tile_url <- function(t) img$getDownloadURL(list(
    crs = "EPSG:4326", format = "GEO_TIFF",
    crs_transform = list(dpp, 0, x0 + t[1] * dpp, 0, -dpp, y0 - t[2] * dpp),
    dimensions = sprintf("%dx%d", t[3], t[4])))
  tile_path <- function(t) file.path(tile_dir, sprintf("t_%d_%d_%d_%d.tif", t[1], t[2], t[3], t[4]))

  # Tiles download in a pool of worker processes that keeps CONCURRENT requests
  # in flight, so one slow tile never holds up the rest. URLs are made here,
  # where the Earth Engine session lives; workers only download them.
  cl <- NULL
  on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
  worker <- fetch_tile
  environment(worker) <- baseenv()
  job_for <- function(t) {
    url <- tryCatch(tile_url(t), error = function(e) e)
    list(url = if (!inherits(url, "error")) url, path = tile_path(t), tries = TRIES,
         refusal_pattern = GEE_REFUSAL_PATTERN,
         error = if (inherits(url, "error")) conditionMessage(url))
  }

  pending <- cut(0L, 0L, W, H, cap_px)
  sdm_info(sprintf("Downloading %d x %d px as %d tile%s", W, H, length(pending),
                   if (length(pending) == 1L) "" else "s"), indent = 2L)
  got <- character(0); refused_at_min <- NULL
  t_start <- proc.time()[["elapsed"]]; next_report <- 0.1
  while (length(pending) && is.null(refused_at_min)) {
    # Tiles larger than the size last refused are split before being asked for.
    pending <- unlist(lapply(pending, function(t)
      if (max(t[3], t[4]) > cap_px) cut(t[1], t[2], t[3], t[4], cap_px) else list(t)),
      recursive = FALSE)
    # A few rounds of work per worker at a time, so URLs are fresh when used and a
    # refusal shrinks the tiles still waiting.
    batch   <- pending[seq_len(min(4L * CONCURRENT, length(pending)))]
    pending <- pending[-seq_along(batch)]
    jobs    <- lapply(batch, job_for)
    results <- if (length(jobs) == 1L) list(worker(jobs[[1]])) else {
      if (is.null(cl)) cl <- parallel::makeCluster(min(CONCURRENT, length(jobs) + length(pending)))
      parallel::parLapplyLB(cl, jobs, worker)
    }
    for (i in seq_along(batch)) {
      t <- batch[[i]]; res <- results[[i]]
      if (isTRUE(res)) { got <- c(got, jobs[[i]]$path); next }
      if (!ee_refused(res))
        stop(sprintf("Map tile download failed after %d tries: %s", TRIES, res), call. = FALSE)
      if (max(t[3], t[4]) <= MIN_TILE_PX) { refused_at_min <- res; break }
      cap_px <- max(MIN_TILE_PX, max(t[3], t[4]) %/% 2L)
      sdm_info(sprintf("A tile was too expensive to compute; continuing at %d px.", cap_px),
               indent = 2L)
      pending <- c(pending, cut(t[1], t[2], t[3], t[4], cap_px))
    }
    frac <- length(got) / (length(got) + length(pending))
    if (length(pending) && frac >= next_report) {
      sdm_info(sprintf("%d tiles done, %d to go (%.0f s)", length(got), length(pending),
                       proc.time()[["elapsed"]] - t_start), indent = 2L)
      next_report <- frac + 0.1
    }
  }

  if (!is.null(refused_at_min)) {
    sdm_info(paste("Earth Engine will not compute this map tile by tile, so it goes to",
                   "the batch system, which writes tiles to Google Drive; AlphaSDM",
                   "downloads and then removes them. This is slower."), indent = 2L)
    unlink(got)
    got <- list.files(ee_export_image_drive(image$toFloat(), region, scale, tile_dir,
                                            nodata = nodata),
                      pattern = "\\.tif$", full.names = TRUE)
    if (!length(got)) stop("The batch export produced no tiles.", call. = FALSE)
  }

  # Stitch through a virtual mosaic, one band at a time, so nothing has to fit
  # in memory. Earth Engine's multi-band GeoTIFFs make libtiff note, for every
  # tile it opens, that non-colour bands are treated as extra samples; that is
  # how the bands should be read, so only that note is silenced.
  quiet_gdal <- function(expr) withCallingHandlers(expr, warning = function(w)
    if (grepl("ExtraSamples", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning"))
  vrt <- file.path(tile_dir, "mosaic.vrt")
  quiet_gdal(sf::gdal_utils("buildvrt", got, vrt,
                            options = c("-srcnodata", nodata, "-vrtnodata", nodata)))
  for (b in seq_along(dsn)) {
    quiet_gdal(sf::gdal_utils("translate", vrt, dsn[b], options = c(
      "-b", b, "-a_nodata", nodata, "-co", "COMPRESS=DEFLATE", "-co", "TILED=YES",
      "-co", "BIGTIFF=IF_SAFER")))
  }
  dsn
}

#' Resolve a user AOI into an ee.Geometry
#'
#' Accepts a pre-built ee.Geometry, a list(lon, lat, radius) buffer, a
#' path to a vector file, or "bbox" for the bounding box of `data`'s points.
#' The file path converts through sf alone (never rgee::sf_as_ee, whose
#' geojsonio dependency is not declared anywhere).
#' @noRd
resolve_aoi <- function(aoi, ee = reticulate::import("ee"), data = NULL) {
  if (inherits(aoi, "python.builtin.object")) return(aoi)
  if (identical(aoi, "bbox") && !is.null(data))
    return(ee$Geometry$Rectangle(c(min(data$longitude), min(data$latitude),
                                   max(data$longitude), max(data$latitude))))
  if (is.list(aoi) && !is.null(aoi$lat))
    return(ee$Geometry$Point(c(as.numeric(aoi$lon), as.numeric(aoi$lat)))$
             buffer(as.numeric(aoi$radius)))
  if (is.character(aoi) && file.exists(aoi)) {
    aoi_sf <- sf::st_read(aoi, quiet = TRUE)
    geom <- sf::st_geometry(aoi_sf)
    if (length(geom) > 1L) geom <- sf::st_union(geom)
    geom <- sf::st_transform(geom, 4326)
    tmp <- tempfile(fileext = ".geojson")
    on.exit(unlink(tmp), add = TRUE)
    sf::st_write(sf::st_sf(geometry = geom), tmp, quiet = TRUE)
    gj <- jsonlite::fromJSON(readLines(tmp, warn = FALSE),
                             simplifyVector = FALSE)
    return(ee$Geometry(gj$features[[1]]$geometry))
  }
  stop("`aoi` must be an ee.Geometry, a list with lon/lat/radius, a path to ",
       "a readable vector file, or \"bbox\". Got: ",
       paste(class(aoi), collapse = "/"), call. = FALSE)
}
