#' Internal Unified GEE Training Pipeline
#' @keywords internal
fit_gee_models <- function(train_df, methods, aoi_geom, scale, aoi_year, training_params, count = NULL,
                           bg_ratio = NULL, persist_classifier = FALSE, project = NULL) {
  ee <- reticulate::import("ee")

  # Optional class balancing: downsample an absence/background FeatureCollection to
  # a target absence:presence ratio, entirely server-side (no embeddings egressed).
  # On imbalanced SDM data this is the main lever that moves TSS for the tree
  # methods (rf/gbt). Only downsamples — never upsamples — so a ratio looser than
  # the data already is leaves the pool untouched.
  balance_bg <- function(bg_fc, n_pos, n_neg) {
    if (is.null(bg_ratio) || n_pos <= 0L || n_neg <= 0L) return(bg_fc)
    target <- ceiling(as.numeric(bg_ratio) * n_pos)
    if (target >= n_neg) return(bg_fc)
    frac <- target / n_neg
    bg_fc$randomColumn("__bgsel", 42L)$filter(ee$Filter$lt("__bgsel", frac))
  }

  # 1. Upload + background generation (background stays on GEE if presence-only)
  sdm_section("Uploading training data to Google Earth Engine")
  pb_up <- sdm_progress_start("Uploading and sampling")

  if (all(train_df$present == 1)) {
    pres_df  <- train_df[, c("longitude", "latitude", "year", "present")]
    n_pres   <- nrow(pres_df)
    bg_count <- if (is.null(count)) n_pres else as.integer(count)

    lons <- pres_df$longitude; lats <- pres_df$latitude
    sdm_info(sprintf("Coordinate extent — Lon [%.2f, %.2f]  Lat [%.2f, %.2f]",
                     min(lons), max(lons), min(lats), max(lats)), indent = 1L)
    sdm_info(sprintf("Transferring %d presence coordinates ...", n_pres), indent = 1L)

    all_years   <- as.list(unique(c(as.integer(pres_df$year), as.integer(aoi_year))))
    presence_fc <- upload_points_to_gee(pres_df)

    pres_sampled <- get_embeddings_at_fc(presence_fc, scale,
                                         properties = c("year", "present"),
                                         years      = all_years)

    # Background: server-side randomPoints — computation graph stays small regardless
    # of n, avoiding the 10 MB payload limit from embedding re-upload.
    lon_buf   <- max(0.1, (max(lons) - min(lons)) * 0.1)
    lat_buf   <- max(0.1, (max(lats) - min(lats)) * 0.1)
    bbox_flat <- c(min(lons)-lon_buf, min(lats)-lat_buf, max(lons)+lon_buf, max(lats)+lat_buf)

    n_bg <- min(n_pres, bg_count)
    sdm_info(sprintf("Sampling %d background points (server-side) ...", n_bg), indent = 1L)
    bg_fc <- generate_background_fc_gee(bbox_flat, aoi_year, n_bg, aoi_geom = aoi_geom)

    bg_props   <- c("year", "present")
    bg_sampled <- get_embeddings_at_fc(bg_fc, scale,
                                       properties = bg_props, years = list(as.integer(aoi_year)))

    bg_rf      <- balance_bg(bg_sampled, n_pres, n_bg)
    bg_gbt     <- bg_rf
    sampled_fc <- pres_sampled$merge(bg_sampled)
    n_background <- n_bg
  } else {
    upload_df    <- train_df[, c("longitude", "latitude", "year", "present")]
    sdm_info(sprintf("Transferring %d coordinates ...", nrow(upload_df)), indent = 1L)

    upload_fc    <- upload_points_to_gee(upload_df)
    sampled_fc   <- get_embeddings_at_fc(upload_fc, scale,
                                         properties = c("year", "present"),
                                         years      = as.list(unique(as.integer(upload_df$year))))
    pres_sampled <- sampled_fc$filter(ee$Filter$eq("present", 1L))
    # When real absences are supplied, RF/GBT use them as background (merged back
    # with presences below) — mirrors the presence-only branch's bg_rf/bg_gbt.
    n_pres       <- sum(train_df$present == 1)
    n_background <- sum(train_df$present == 0)
    bg_rf        <- balance_bg(sampled_fc$filter(ee$Filter$eq("present", 0L)), n_pres, n_background)
    bg_gbt       <- bg_rf
  }

  sdm_progress_done(pb_up)

  # 2. Train all models, forcing each classifier to evaluate immediately (chunked GEE training).
  # GEE classifiers are lazy: clf$train() only builds a computation graph.
  # By calling getInfo() on each classifier right after training, we materialize them
  # one at a time — one getInfo() per model keeps each call within GEE's timeout.
  # Similarity is already eager (reduceColumns$getInfo() fires inside train_gee_model).
  sdm_section(sprintf("Training %d model%s on Google Earth Engine",
                      length(methods), if (length(methods) == 1) "" else "s"))
  pb <- sdm_progress_start("Model training", total = length(methods))
  models <- list()
  for (m in methods) {
    sdm_info(sprintf("Fitting %s ...", toupper(m)), indent = 1L)
    fc_for_method <- if (m == "similarity") {
      pres_sampled
    } else if (m == "rf") {
      pres_sampled$merge(bg_rf)
    } else if (m == "gbt") {
      pres_sampled$merge(bg_gbt)
    } else {
      sampled_fc
    }
    models[[m]]   <- train_gee_model(fc_for_method, m, params = training_params[[m]],
                                     persist = persist_classifier, project = project)

    if (models[[m]]$is_classifier) {
      invisible(retry_curl_download(
        pres_sampled$limit(1L)$classify(models[[m]]$trained)$getInfo()
      ))
    }

    pb <- sdm_progress_update(pb)
  }
  sdm_progress_done(pb)

  return(list(
    models   = models,
    metadata = list(
      n_presence       = n_pres,
      n_background     = n_background,
      methods          = methods,
      scale            = scale,
      # temporary classifier assets created by persist_classifier (NULL if none),
      # to be removed with cleanup_classifier_assets() once prediction is done.
      classifier_assets = Filter(Negate(is.null), lapply(models, function(x) x$asset_id))
    )
  ))
}

#' Delete temporary classifier assets created by persist_classifier
#' @keywords internal
cleanup_classifier_assets <- function(train_res) {
  assets <- train_res$metadata$classifier_assets
  if (length(assets) > 0) for (a in assets) ee_delete_asset_quietly(a)
  invisible(NULL)
}


#' Generate SDM Map
#' @export
generate_map <- function(data, aoi, scale = 10, output_dir = getwd(),
                         methods = NULL, ensemble = TRUE, aoi_year = NULL, count = NULL, bg_ratio = NULL,
                         n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                         shrinkage = 0.005, max_nodes = 6L, variables_per_split = NULL,
                         svm_type = "EPSILON_SVR", svm_kernel = "RBF", svm_cost = 10, svm_gamma = 0.05,
                         maxent_beta = 1, maxent_features = "auto",
                         persist_classifier = FALSE,
                         gee_project = NULL, python_path = NULL,
                         options = list()) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)
  t_total_start <- proc.time()[["elapsed"]]

  ee <- reticulate::import("ee")

  if (is.null(aoi_year)) aoi_year <- 2023
  # Default ensemble: the strong, complementary tier validated across the
  # benchmarks. Lighter models (similarity, knn, cart, mindist) remain available
  # via `methods=` but trail by ~0.03-0.07 AUC, so they are not in the default.
  if (is.null(methods)) methods <- c("svm", "rf", "gbt", "maxent")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. Prepare AOI Geometry
  aoi_geom <- NULL
  if (inherits(aoi, "python.builtin.object")) {
    aoi_geom <- aoi                                   # pre-built ee.Geometry
  } else if (is.list(aoi) && !is.null(aoi$lat)) {
    aoi_geom <- ee$Geometry$Point(c(as.numeric(aoi$lon), as.numeric(aoi$lat)))$buffer(as.numeric(aoi$radius))
  } else if (is.character(aoi) && file.exists(aoi)) {
    aoi_sf <- sf::st_read(aoi, quiet = TRUE)
    aoi_geom <- rgee::sf_as_ee(aoi_sf)$geometry()
  }

  # 2. Method-Specific Optimized Defaults
  base_params <- list(
    numberOfTrees = n_trees, 
    minLeafPopulation = min_leaf_population,
    bagFraction = bag_fraction, 
    shrinkage = shrinkage,
    maxNodes = max_nodes, 
    variablesPerSplit = variables_per_split
  )

  # Setup method-specific params list
  method_params <- list()
  for (m in methods) {
    method_params[[m]] <- base_params
  }

  # Override GBT with 150 trees if using global default of 100
  if ("gbt" %in% methods && n_trees == 100L) {
    method_params$gbt$numberOfTrees <- 150L
  }
  
  # Override RF with 250 trees if using global default of 100
  if ("rf" %in% methods && n_trees == 100L) {
    method_params$rf$numberOfTrees <- 250L
  }
  # RF benefits from DEEP trees; the shared default maxNodes=6 / minLeaf=5 are near-stumps
  # (benchmarked +0.02-0.05 AUC from unlimited depth). RF-specific so gbt/cart stay shallow.
  if ("rf" %in% methods && max_nodes == 6L)            method_params$rf$maxNodes <- NULL
  if ("rf" %in% methods && min_leaf_population == 5L)  method_params$rf$minLeafPopulation <- 1L

  # SVM tuning: benchmarked ε-SVR + RBF defaults (overridable per call)
  if ("svm" %in% methods) {
    method_params$svm$svmType    <- svm_type
    method_params$svm$kernelType <- svm_kernel
    method_params$svm$cost       <- svm_cost
    method_params$svm$gamma      <- svm_gamma
  }
  # MaxEnt tuning (ENMeval-style): regularization multiplier + feature classes
  if ("maxent" %in% methods)
    method_params$maxent <- modifyList(method_params$maxent, maxent_tuning_params(maxent_beta, maxent_features))

  # 4. Unified Training
  train_res <- fit_gee_models(data, methods, aoi_geom, scale, aoi_year, method_params, count = count, bg_ratio = bg_ratio, persist_classifier = persist_classifier, project = gee_project)
  on.exit(cleanup_classifier_assets(train_res), add = TRUE)   # remove temp classifier assets on exit

  # 5. Map Generation
  img_mosaic <- get_embedding_image(aoi_year, scale)
  final_results <- list(methods = methods, model_metadata = train_res$metadata)
  pb_map <- sdm_progress_start("Map generation", total = length(methods))
  for (m in methods) {
    sdm_info(sprintf("Exporting %s ...", toupper(m)), indent = 1L)
    pred_img <- predict_gee_map(train_res$models[[m]], img_mosaic)
    tif_path <- file.path(output_dir, paste0(m, ".tif"))
    # Drive-free, water-robust export: tile getDownloadURL + local mosaic.
    tryCatch(
      export_image_tiled(pred_img, aoi_geom, scale, tif_path),
      error = function(e) sdm_warn(sprintf("Export of %s failed: %s", m, conditionMessage(e)))
    )
    final_results[[paste0(m, "_map")]] <- tif_path
    pb_map <- sdm_progress_update(pb_map)
  }
  sdm_progress_done(pb_map)
  
  t_total_end <- proc.time()[["elapsed"]]
  sdm_section("Species processing finished")
  sdm_done(sprintf("Total elapsed time [%.1fs]", t_total_end - t_total_start))
  cat("\n")
  
  # Reset standardization trackers for next run
  .alphasdm_env$standardization_active <- FALSE
  .alphasdm_env$standardization_info_printed <- FALSE
  
  return(final_results)
}

#' Internal: score a coordinate data frame against trained models (image-first, server-side)
#'
#' Classifies the embedding image once per model, stacks into a multi-band prediction image,
#' then batch-samples that image at eval coords. This avoids re-embedding eval points through
#' 64 bands and re-triggering classifier training on every batch.
#' @keywords internal
predict_scores_internal <- function(predict_df, models, methods, img, scale, aoi_year,
                                    async = FALSE, project = NULL) {
  ee <- reticulate::import("ee")
  if (!"year" %in% names(predict_df)) predict_df$year <- aoi_year

  # Build multi-band prediction image once — classifier training evaluated here by GEE,
  # then cached for all subsequent sampleRegions batches.
  pred_imgs  <- lapply(methods, function(m) predict_gee_map(models[[m]], img)$rename(paste0("pred_", m)))
  pred_stack <- Reduce(function(a, b) a$addBands(b), pred_imgs)

  pred_cols <- paste0("pred_", methods)
  for (col in pred_cols) predict_df[[col]] <- NA_real_
  predict_df$.row_idx <- seq_len(nrow(predict_df))

  fill_from_features <- function(features) {
    for (f in features) {
      ridx <- as.integer(f$properties[["row_idx"]])
      if (is.null(ridx) || is.na(ridx)) next
      for (m in methods) {
        val <- f$properties[[paste0("pred_", m)]]
        if (!is.null(val)) predict_df[ridx, paste0("pred_", m)] <<- as.numeric(val)
      }
    }
  }

  # --- Async path: one batch export of the whole scored collection -------------
  # Materialises training + sampling server-side (no synchronous compute limit),
  # then reads the result back in pages. Embeddings never leave GEE.
  if (async) {
    chunk <- predict_df[, c("longitude", "latitude", "year", ".row_idx"), drop = FALSE]
    names(chunk)[names(chunk) == ".row_idx"] <- "row_idx"
    eval_fc    <- upload_points_to_gee(chunk)
    # geometries = TRUE: Export.table.toAsset rejects null-geometry features.
    sampled_fc <- pred_stack$sampleRegions(
      collection = eval_fc, properties = as.list("row_idx"),
      scale = as.integer(scale), tileScale = 16L, geometries = TRUE
    )
    info <- ee_table_to_info_async(sampled_fc$select(as.list(c("row_idx", pred_cols))), project)
    fill_from_features(info$features)
    predict_df$.row_idx <- NULL
    return(predict_df)
  }

  # --- Synchronous path (default): chunked getInfo -----------------------------
  chunk_size   <- 5000L
  n_rows       <- nrow(predict_df)
  batch_starts <- seq(1L, n_rows, by = chunk_size)

  for (i in batch_starts) {
    end_idx <- min(i + chunk_size - 1L, n_rows)
    chunk   <- predict_df[i:end_idx, c("longitude", "latitude", "year", ".row_idx"), drop = FALSE]
    names(chunk)[names(chunk) == ".row_idx"] <- "row_idx"

    eval_fc    <- upload_points_to_gee(chunk)
    sampled_fc <- pred_stack$sampleRegions(
      collection = eval_fc,
      properties = as.list("row_idx"),
      scale      = as.integer(scale),
      tileScale  = 16L,
      geometries = FALSE
    )

    res <- tryCatch(retry_curl_download(sampled_fc$getInfo()), error = function(e) e)
    if (inherits(res, "error")) {
      # A genuine compute timeout means the job is too big for synchronous scoring —
      # fail loudly and point the user at the async path rather than silently
      # returning NA scores (which would corrupt downstream metrics).
      if (is_gee_timeout(res)) stop(async_timeout_message(res, where = "scoring"))
      next  # other transient errors: skip this chunk (legacy behaviour)
    }
    if (is.null(res) || !"features" %in% names(res)) next
    fill_from_features(res$features)
  }

  predict_df$.row_idx <- NULL
  return(predict_df)
}


#' Evaluate SDM Models
#' @export
evaluate_models <- function(data, predict_coords = NULL, scale = 10, output_dir = getwd(),
                            methods = NULL, aoi_year = NULL, count = NULL, bg_ratio = NULL,
                            n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                            shrinkage = 0.005, max_nodes = 6L, variables_per_split = NULL,
                            svm_type = "EPSILON_SVR", svm_kernel = "RBF", svm_cost = 10, svm_gamma = 0.05,
                            maxent_beta = 1, maxent_features = "auto",
                            cv_folds = 5L, cv_method = "block", block_size = NULL,
                            weighted_ensemble = FALSE, async = FALSE,
                            persist_classifier = FALSE,
                            gee_project = NULL, python_path = NULL,
                            options = list()) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)
  t_total_start <- proc.time()[["elapsed"]]

  ee <- reticulate::import("ee")

  if (is.null(aoi_year)) aoi_year <- 2023
  # Default ensemble: the strong, complementary tier validated across the
  # benchmarks. Lighter models (similarity, knn, cart, mindist) remain available
  # via `methods=` but trail by ~0.03-0.07 AUC, so they are not in the default.
  if (is.null(methods)) methods <- c("svm", "rf", "gbt", "maxent")

  # --- Parameter validation ---
  cv_folds <- as.integer(cv_folds)
  if (weighted_ensemble && cv_folds == 0L) {
    cv_folds <- 5L
    sdm_warn("weighted_ensemble = TRUE requires cross-validation; cv_folds overridden to 5.")
  }
  if (is.null(predict_coords) && cv_folds == 0L) {
    stop("Either predict_coords or cv_folds > 0 must be provided.")
  }

  # --- Method-Specific Optimized Defaults ---
  base_params <- list(
    numberOfTrees = n_trees, minLeafPopulation = min_leaf_population,
    bagFraction = bag_fraction, shrinkage = shrinkage,
    maxNodes = max_nodes, variablesPerSplit = variables_per_split
  )
  method_params <- setNames(lapply(methods, function(m) base_params), methods)
  # Validated defaults: gbt 150 trees @ shrinkage 0.05, maxNodes 6;
  # rf 500 deep trees, variablesPerSplit 8. Each override only fires when the corresponding
  # global argument is still at its default, so explicit user values always win.
  if ("gbt" %in% methods && n_trees == 100L)    method_params$gbt$numberOfTrees <- 150L
  if ("gbt" %in% methods && shrinkage == 0.005) method_params$gbt$shrinkage     <- 0.05
  if ("rf"  %in% methods && n_trees == 100L)    method_params$rf$numberOfTrees  <- 500L
  if ("rf"  %in% methods && is.null(variables_per_split)) method_params$rf$variablesPerSplit <- 8L
  # RF benefits from deep trees (default maxNodes=6/minLeaf=5 are near-stumps); rf-specific.
  if ("rf"  %in% methods && max_nodes == 6L)           method_params$rf$maxNodes <- NULL
  if ("rf"  %in% methods && min_leaf_population == 5L) method_params$rf$minLeafPopulation <- 1L
  if ("svm" %in% methods) {
    method_params$svm$svmType    <- svm_type
    method_params$svm$kernelType <- svm_kernel
    method_params$svm$cost       <- svm_cost
    method_params$svm$gamma      <- svm_gamma
  }
  # MaxEnt tuning (ENMeval-style): regularization multiplier + feature classes
  if ("maxent" %in% methods)
    method_params$maxent <- modifyList(method_params$maxent, maxent_tuning_params(maxent_beta, maxent_features))

  # --- AOI geometry (bounding box of prediction target or training data) ---
  ref_df   <- if (!is.null(predict_coords)) predict_coords else data
  bbox     <- c(min(ref_df$longitude), min(ref_df$latitude),
                max(ref_df$longitude), max(ref_df$latitude))
  aoi_geom <- ee$Geometry$Rectangle(bbox)

  # Embedding image — lazy GEE object, shared across training and prediction.
  img <- get_embedding_image(aoi_year, scale)

  # --- Cross-validation ---
  cv_metrics       <- NULL
  ensemble_weights <- NULL

  if (cv_folds > 0L) {
    sdm_section(sprintf("Cross-validation (%d spatial folds)", cv_folds))

    pres_df <- data[data$present == 1, , drop = FALSE]
    n_pres  <- nrow(pres_df)

    if (cv_folds > n_pres) {
      stop(sprintf(
        "cv_folds (%d) exceeds the number of presence points (%d). Reduce cv_folds or provide more training data.",
        cv_folds, n_pres
      ))
    }

    # Spatial folds: blockCV blocks by default (controls spatial autocorrelation
    # better than contiguous k-means clusters); falls back to k-means if unavailable.
    pres_df$cv_fold <- assign_cv_folds(pres_df, cv_folds, method = cv_method, block_size = block_size)
    sdm_info(sprintf("CV folds: %s (%d folds)",
                     switch(cv_method, block = "blockCV blocks", kmeans = "k-means clusters",
                            random = "random k-fold", cv_method), cv_folds), indent = 1L)

    # Generate validation background on GEE, download coordinates once before the CV loop
    val_bg_n  <- if (is.null(count)) n_pres else as.integer(count)
    bg_fc_gee <- generate_background_fc_gee(bbox, aoi_year, val_bg_n)
    bg_info   <- retry_curl_download(bg_fc_gee$getInfo())
    val_bg_df <- do.call(rbind, lapply(bg_info$features, function(f) {
      data.frame(longitude = f$geometry$coordinates[[1]],
                 latitude  = f$geometry$coordinates[[2]],
                 year      = aoi_year)
    }))

    cv_fold_aucs <- matrix(NA_real_, nrow = cv_folds, ncol = length(methods),
                           dimnames = list(NULL, methods))

    for (k in seq_len(cv_folds)) {
      sdm_info(sprintf("Fold %d / %d", k, cv_folds), indent = 1L)

      train_k <- pres_df[pres_df$cv_fold != k, c("longitude", "latitude", "year", "present"), drop = FALSE]
      test_k  <- pres_df[pres_df$cv_fold == k, c("longitude", "latitude", "year"),             drop = FALSE]

      res_k        <- fit_gee_models(train_k, methods, aoi_geom, scale, aoi_year, method_params, count = count, bg_ratio = bg_ratio, persist_classifier = persist_classifier, project = gee_project)
      pres_scored  <- predict_scores_internal(test_k,    res_k$models, methods, img, scale, aoi_year, async = async, project = gee_project)
      bg_scored    <- predict_scores_internal(val_bg_df, res_k$models, methods, img, scale, aoi_year, async = async, project = gee_project)
      cleanup_classifier_assets(res_k)   # per-fold temp classifier assets

      na_pres_cv <- is.na(pres_scored[[paste0("pred_", methods[1])]])
      na_bg_cv   <- is.na(bg_scored[[paste0("pred_",  methods[1])]])
      if (any(na_pres_cv) || any(na_bg_cv)) {
        sdm_warn(sprintf(
          "CV fold %d: %d presence and %d background points dropped — no satellite coverage.",
          k, sum(na_pres_cv), sum(na_bg_cv)
        ), indent = 2L)
      }
      for (m in methods) {
        pos <- pres_scored[[paste0("pred_", m)]]; pos <- pos[!is.na(pos)]
        neg <- bg_scored[[paste0("pred_", m)]];   neg <- neg[!is.na(neg)]
        if (length(pos) > 0 && length(neg) > 0)
          cv_fold_aucs[k, m] <- calculate_classifier_metrics(pos, neg)$auc_roc
      }
    }

    cv_aucs    <- colMeans(cv_fold_aucs, na.rm = TRUE)
    cv_metrics <- as.list(cv_aucs)
    sdm_info(paste("CV AUC —", paste(sprintf("%s:%.3f", names(cv_aucs), cv_aucs), collapse = "  ")))

    if (weighted_ensemble) {
      raw_wts <- pmax(cv_aucs - 0.5, 0)
      if (sum(raw_wts) == 0) {
        sdm_warn("All CV AUCs <= 0.5; falling back to equal weights")
        ensemble_weights <- setNames(rep(1 / length(methods), length(methods)), methods)
      } else {
        ensemble_weights <- raw_wts / sum(raw_wts)
      }
      sdm_info(paste("Ensemble weights —", paste(sprintf("%s:%.3f", names(ensemble_weights), ensemble_weights), collapse = "  ")))
    }
  }

  # --- Final model trained on all data ---
  train_res <- fit_gee_models(data, methods, aoi_geom, scale, aoi_year, method_params, count = count, bg_ratio = bg_ratio, persist_classifier = persist_classifier, project = gee_project)
  on.exit(cleanup_classifier_assets(train_res), add = TRUE)   # remove temp assets on any exit

  # --- Pure CV mode: return without predict_coords ---
  if (is.null(predict_coords)) {
    t_total_end <- proc.time()[["elapsed"]]
    sdm_section("Species evaluation finished")
    sdm_done(sprintf("Total elapsed time [%.1fs]", t_total_end - t_total_start))
    cat("\n")
    .alphasdm_env$standardization_active <- FALSE
    .alphasdm_env$standardization_info_printed <- FALSE
    return(list(
      methods          = c(methods, "ensemble"),
      model_metadata   = train_res$metadata,
      cv_metrics       = cv_metrics,
      ensemble_weights = ensemble_weights
    ))
  }

  # --- Prediction ---
  sdm_section(sprintf("Predicting at %d coordinates (image-first, server-side)", nrow(predict_coords)))
  pb_pred    <- sdm_progress_start("Prediction")
  final_pred <- predict_scores_internal(predict_coords, train_res$models, methods, img, scale, aoi_year, async = async, project = gee_project)
  sdm_progress_done(pb_pred)

  # --- Ensemble ---
  pred_cols <- paste0("pred_", methods)
  if (weighted_ensemble && !is.null(ensemble_weights)) {
    wts <- ensemble_weights[methods]
    final_pred$pred_ensemble <- as.numeric(as.matrix(final_pred[, pred_cols, drop = FALSE]) %*% wts)
  } else {
    final_pred$pred_ensemble <- rowMeans(final_pred[, pred_cols, drop = FALSE], na.rm = TRUE)
  }

  # --- Metrics ---
  final_results <- list(
    methods           = c(methods, "ensemble"),
    metrics           = list(),
    point_predictions = final_pred,
    model_metadata    = train_res$metadata,
    cv_metrics        = cv_metrics,
    ensemble_weights  = ensemble_weights
  )

  if ("present" %in% names(final_pred)) {
    # Warn about eval points dropped due to missing satellite coverage.
    # All models share the same image mask, so NA in any one model means NA in all.
    na_mask <- is.na(final_pred[[paste0("pred_", methods[1])]])
    if (any(na_mask)) {
      na_pres <- sum(na_mask & final_pred$present == 1)
      na_bg   <- sum(na_mask & final_pred$present == 0)
      sdm_warn(sprintf(
        "%d eval point%s dropped — no satellite coverage at prediction time (%d presence, %d background).",
        sum(na_mask), if (sum(na_mask) == 1) "" else "s", na_pres, na_bg
      ), indent = 1L)
      dropped_coords <- final_pred[na_mask, c("longitude", "latitude", "present"), drop = FALSE]
      dropped_coords$class <- ifelse(dropped_coords$present == 1, "presence", "background")
      for (j in seq_len(nrow(dropped_coords))) {
        sdm_info(sprintf("  Dropped %s: lon=%.4f, lat=%.4f",
                         dropped_coords$class[j],
                         dropped_coords$longitude[j],
                         dropped_coords$latitude[j]), indent = 2L)
      }
    }

    for (m in c(methods, "ensemble")) {
      scores <- final_pred[[paste0("pred_", m)]]
      pos <- scores[final_pred$present == 1]
      neg <- scores[final_pred$present == 0]
      final_results$metrics[[m]] <- calculate_classifier_metrics(pos, neg)
    }
  }

  t_total_end <- proc.time()[["elapsed"]]
  sdm_section("Species evaluation finished")
  sdm_done(sprintf("Total elapsed time [%.1fs]", t_total_end - t_total_start))
  cat("\n")

  .alphasdm_env$standardization_active <- FALSE
  .alphasdm_env$standardization_info_printed <- FALSE

  return(final_results)
}
