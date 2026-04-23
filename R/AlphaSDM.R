#' Generate SDM Map
#'
#' Trains models on provided data and generates spatial prediction maps (rasters) on GEE.
#'
#' @param data A data frame formatted via `format_data()`. Must have (longitude, latitude, year, present).
#'   If a 'species' column is present with multiple species, the function processes each separately.
#' @param aoi Either a list with `lat`, `lon`, and `radius` (in meters), or a path to a polygon file.
#' @param scale Resolution in meters for the final map. Defaults to 10.
#' @param output_dir Directory to save results. Defaults to the current working directory.
#' @param methods Character vector of method names. Supported: "rf", "gbt", "maxent", "svm", "ridge", "centroid".
#' @param ensemble Combine all methods into an ensemble map. Defaults to TRUE.
#' @param aoi_year Alpha Earth Mosaic year for map generation. Defaults to current year.
#' @param count Number of background points. Defaults to 10x presence points.
#' @param n_trees Number of trees for tree-based methods.
#' @param min_leaf_population Min leaf size for trees.
#' @param bag_fraction Bagging fraction.
#' @param shrinkage Learning rate for GBT.
#' @param max_nodes Max nodes per tree.
#' @param variables_per_split Variables per split.
#' @param lambda_ Penalty for Ridge.
#' @param polynomial Degree of polynomial terms for Ridge regression. 1 = linear, 2 = quadratic (default).
#' @param gee_project Optional override. Normally configured once via \code{\link{setup_gee}}.
#' @param python_path Optional. Path to Python executable.
#' @param options Named list of advanced hyperparameters.
#'
#' @importFrom stats setNames runif complete.cases na.omit
#' @importFrom utils read.csv
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply detectCores
#' @export
generate_map <- function(data, aoi, scale = 10, output_dir = getwd(),
                         methods = NULL, ensemble = TRUE, aoi_year = NULL, count = NULL,
                         n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                         shrinkage = 0.005, max_nodes = NULL, variables_per_split = NULL, lambda_ = 0.0,
                         polynomial = 2L,
                         gee_project = NULL, python_path = NULL,
                         options = list()) {

    if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)
  ee <- reticulate::import("ee")

  # Multi-species check
  if ("species" %in% names(data) && length(unique(data$species)) > 1) {
    species_list <- unique(data$species)
    timestamp_message(sprintf("--- Multi-Species Mode: Processing %d species ---", length(species_list)))
    results_list <- list()
    for (sp in species_list) {
      timestamp_message(sprintf("\n--- Species: %s ---", sp))
      results_list[[sp]] <- generate_map(
        data[data$species == sp, ], aoi, scale, file.path(output_dir, sp),
        methods, ensemble, aoi_year, count, n_trees, min_leaf_population,
        bag_fraction, shrinkage, max_nodes, variables_per_split, lambda_,
        polynomial,
        gee_project, python_path, options
      )
    }
    return(results_list)
  }

  if (is.null(aoi_year)) aoi_year <- 2023
  if (is.null(methods)) methods <- c("centroid", "ridge", "rf", "gbt", "maxent", "svm")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. Prepare AOI Geometry
  aoi_geom <- NULL
  if (is.list(aoi) && !is.null(aoi$lat)) {
    aoi_geom <- ee$Geometry$Point(c(as.numeric(aoi$lon), as.numeric(aoi$lat)))$buffer(as.numeric(aoi$radius))
  } else if (is.character(aoi) && file.exists(aoi)) {
    aoi_sf <- sf::st_read(aoi, quiet = TRUE)
    aoi_geom <- rgee::sf_as_ee(aoi_sf)$geometry()
  } else {
    stop("Invalid AOI. Please provide a center/radius list or a valid spatial file path.")
  }

  # 2. Train Models (Shared Pipeline)
  training_params <- list(
    numberOfTrees = n_trees, minLeafPopulation = min_leaf_population,
    bagFraction = bag_fraction, shrinkage = shrinkage,
    maxNodes = max_nodes, variablesPerSplit = variables_per_split, lambda_ = lambda_,
    polynomial = polynomial
  )
  for (opt in names(options)) training_params[[opt]] <- options[[opt]]

  train_res <- train_AlphaSDM_internal(
    data = data, methods = methods, aoi_geom = aoi_geom, scale = scale,
    aoi_year = aoi_year, count = count, training_params = training_params
  )

  models <- train_res$models
  final_results <- list(methods = methods, model_metadata = train_res$metadata)

  # 3. Map Generation
  img_mosaic <- get_embedding_image(aoi_year, scale)

  # Individual maps
  timestamp_message("--- Generating Maps on GEE ---")
  for (m in methods) {
    pred_img <- predict_gee_map(models[[m]], img_mosaic)
    tif_path <- normalizePath(file.path(output_dir, paste0(gsub("[^a-zA-Z0-9]", "", m), ".tif")), mustWork = FALSE)
    timestamp_message(sprintf("  Saving map to: %s", tif_path))
    try(rgee::ee_as_rast(image = pred_img, region = aoi_geom$bounds(), scale = scale, dsn = tif_path), silent = TRUE)
    final_results[[paste0(m, "_map")]] <- tif_path
  }

  # Ensemble map
  if (ensemble && length(methods) > 1) {
    timestamp_message("  Calculating ensemble map...")
    pred_images <- lapply(models, function(mod) predict_gee_map(mod, img_mosaic))
    ensemble_img <- ee$ImageCollection(unname(pred_images))$mean()$rename("similarity")
    tif_path <- file.path(output_dir, "ensemble.tif")
    try(rgee::ee_as_rast(image = ensemble_img, region = aoi_geom$bounds(), scale = scale, dsn = tif_path), silent = TRUE)
    final_results$ensemble_map <- tif_path
  }

  return(final_results)
}

#' Evaluate SDM Models
#'
#' @param data A data frame formatted via `format_data()`.
#' @param predict_coords Data frame of target coordinates.
#' @param methods Character vector of method names.
#' @param scale Resolution in meters.
#' @param output_dir Directory to save results.
#' @param aoi_year Alpha Earth Mosaic year.
#' @param count Number of background points.
#' @param n_trees Number of trees for tree-based methods.
#' @param min_leaf_population Min leaf size for trees.
#' @param bag_fraction Bagging fraction.
#' @param shrinkage Learning rate for GBT.
#' @param max_nodes Max nodes per tree.
#' @param variables_per_split Variables per split.
#' @param lambda_ Penalty for Ridge.
#' @param polynomial Degree of polynomial terms for Ridge regression.
#' @param gee_project Optional override. Normally configured once via \code{\link{setup_gee}}.
#' @param python_path Optional. Path to Python executable.
#' @param options Named list of advanced hyperparameters.
#'
#' @export
evaluate_models <- function(data, predict_coords, scale = 10, output_dir = getwd(),
                            methods = NULL, aoi_year = NULL, count = NULL,
                            n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                            shrinkage = 0.005, max_nodes = NULL, variables_per_split = NULL, lambda_ = 0.0,
                            polynomial = 2L,
                            gee_project = NULL, python_path = NULL,
                            options = list()) {

    if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  
  # --- 0. PRE-FLIGHT VALIDATION (Catch naming errors early) ---
  req_cols <- c("longitude", "latitude", "year", "present")
  missing_data <- setdiff(req_cols, names(data))
  if (length(missing_data) > 0) {
      # Help the user catch 'presence' vs 'present'
      if ("presence" %in% names(data)) {
          stop("Naming error: AlphaSDM expects 'present' column, but found 'presence'. Please use 'format_data()' to standardize your inputs.")
      }
      stop(sprintf("Training data is missing required columns: %s", paste(missing_data, collapse = ", ")))
  }
  
  # Validate predict_coords (if provided)
  predict_req <- c("longitude", "latitude", "year")
  missing_predict <- setdiff(predict_req, names(predict_coords))
  if (length(missing_predict) > 0) {
      stop(sprintf("Prediction coordinates are missing required columns: %s", paste(missing_predict, collapse = ", ")))
  }
  
  # Fuzzy check for 'presence' in predict_coords (common failure)
  if ("presence" %in% names(predict_coords) && !"present" %in% names(predict_coords)) {
      stop("Naming error: 'presence' column detected in predict_coords. For AlphaSDM to calculate metrics, the target column must be named 'present'. Use 'format_data()'.")
  }

  ensure_gee_authenticated(project = gee_project)
  ee <- reticulate::import("ee")

  # Multi-species check
  if ("species" %in% names(data) && length(unique(data$species)) > 1) {
    species_list <- unique(data$species)
    timestamp_message(sprintf("--- Multi-Species Mode: Evaluating %d species ---", length(species_list)))
    results_list <- list()
    for (sp in species_list) {
      timestamp_message(sprintf("\n--- Species: %s ---", sp))
      sp_predict <- if ("species" %in% names(predict_coords)) predict_coords[predict_coords$species == sp, ] else predict_coords
      results_list[[sp]] <- evaluate_models(
        data[data$species == sp, ], sp_predict, methods, scale, file.path(output_dir, sp),
        aoi_year, count, n_trees, min_leaf_population,
        bag_fraction, shrinkage, max_nodes, variables_per_split, lambda_,
        polynomial,
        gee_project, python_path, options
      )
    }
    return(results_list)
  }

  if (is.null(aoi_year)) aoi_year <- 2023
  if (is.null(methods)) methods <- c("centroid", "ridge", "rf", "gbt", "maxent", "svm")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. Prepare AOI Geometry for Background Generation (Bounding Box of predict_coords)
  bbox <- c(
    min(predict_coords$longitude), min(predict_coords$latitude),
    max(predict_coords$longitude), max(predict_coords$latitude)
  )
  aoi_geom <- ee$Geometry$Rectangle(bbox)

  # 2. Train Models (Shared Pipeline)
  training_params <- list(
    numberOfTrees = n_trees, minLeafPopulation = min_leaf_population,
    bagFraction = bag_fraction, shrinkage = shrinkage,
    maxNodes = max_nodes, variablesPerSplit = variables_per_split, lambda_ = lambda_,
    polynomial = polynomial
  )
  for (opt in names(options)) training_params[[opt]] <- options[[opt]]

  train_res <- train_AlphaSDM_internal(
    data = data, methods = methods, aoi_geom = aoi_geom, scale = scale,
    aoi_year = aoi_year, count = count, training_params = training_params
  )

  models <- train_res$models
  final_results <- list(methods = methods, model_metadata = train_res$metadata)

  # 3 & 4. Batch Coordinate Prediction & Fast Egress
  timestamp_message(sprintf("--- Predicting at %d coordinates (Server-side, Batched) ---", nrow(predict_coords)))

  selectors <- c("tmp_id", paste0("pred_", methods))
  if (length(methods) > 1) selectors <- c(selectors, "pred_ensemble")

  # Prediction batch size (Default 5000, overridable by options)
  chunk_size <- if (!is.null(options$batch_size)) {
      as.integer(options$batch_size)
  } else {
      max(5000L, ceiling(nrow(predict_coords) / 40))
  }

  batch_starts <- seq(1, nrow(predict_coords), by = chunk_size)
  timestamp_message(sprintf(
    "  Processing %d coordinates across %d parallel batches (Batch size: %d)...",
    nrow(predict_coords), length(batch_starts), chunk_size
  ))

  # Run full prediction pipeline in parallel via mclapply (fork-safe on Linux).
  # Each worker independently: uploads coords → samples embeddings → runs models → downloads CSV.
  # This is the fastest approach because getDownloadURL() triggers the actual GEE computation,
  # so it MUST happen inside the parallel workers to achieve true concurrency.
  n_cores <- min(length(batch_starts), parallel::detectCores(logical = FALSE), 40L)
  timestamp_message(sprintf("  -> Running %d batches in parallel (%d cores)...", length(batch_starts), n_cores))
  
  process_batch <- function(idx) {
    i <- batch_starts[idx]
    end_idx <- min(i + chunk_size - 1, nrow(predict_coords))
    chunk_df <- predict_coords[i:end_idx, , drop = FALSE]

    # 1. Base Sampling (Common for all models)
    upload_fc <- upload_points_to_gee(chunk_df, id_offset = i - 1)
    sampled_fc <- get_embeddings_at_fc(upload_fc, scale, properties = c("year", "tmp_id"))

    # 2. Resilient Model Prediction (Decoupled per method)
    chunk_results <- data.frame(tmp_id = (i - 1):(end_idx - 1))
    
    for (m in methods) {
      m_col <- paste0("pred_", m)
      model_obj <- models[m] # models is a list
      scored_fc <- predict_all_models_gee(sampled_fc, model_obj)
      
      max_retries <- 5L
      model_df <- NULL
      for (attempt in seq_len(max_retries)) {
        res <- tryCatch({
          url <- scored_fc$getDownloadURL(filetype = "CSV", selectors = as.list(c("tmp_id", m_col)))
          read.csv(url)
        }, error = function(e) e)
        
        if (is.data.frame(res)) {
          model_df <- res
          break
        }
        
        if (attempt == max_retries) {
          cat(sprintf("\n  ! Warning: Model %s failed for batch %d: %s\n", m, idx, conditionMessage(res)))
        } else {
          Sys.sleep(2^(attempt - 1))
        }
      }
      
      if (!is.null(model_df)) {
        chunk_results <- merge(chunk_results, model_df, by = "tmp_id", all.x = TRUE)
      } else {
        chunk_results[[m_col]] <- as.numeric(NA)
      }
    }

    # 3. Ensemble calculation (R-side for resilience)
    if (length(methods) > 1) {
      pred_cols <- paste0("pred_", methods)
      # Calculate means excluding NAs where possible
      chunk_results$pred_ensemble <- rowMeans(chunk_results[, pred_cols, drop = FALSE], na.rm = TRUE)
      # If all models failed, result will be NaN. Convert to NA.
      chunk_results$pred_ensemble[is.nan(chunk_results$pred_ensemble)] <- NA
    }

    return(chunk_results)
  }
  
  t_start_pred <- Sys.time()
  if (.Platform$OS.type == "unix") {
    all_res_list <- parallel::mclapply(seq_along(batch_starts), process_batch, mc.cores = n_cores)
  } else {
    all_res_list <- lapply(seq_along(batch_starts), process_batch)
  }
  t_end_pred <- Sys.time()
  timestamp_message(sprintf("  -> All %d batches complete in %.1f seconds.", length(batch_starts), 
                  as.numeric(difftime(t_end_pred, t_start_pred, units = "secs"))))

  res_df <- do.call(rbind, all_res_list)

  # Re-attach scores to predict_coords safely
  final_pred <- predict_coords
  for (col in setdiff(selectors, "tmp_id")) {
    final_pred[[col]] <- NA # initialize natively
  }

  if (!is.null(res_df) && nrow(res_df) > 0) {
    res_df <- res_df[order(res_df$tmp_id), ]
    idx <- match(res_df$tmp_id, (seq_len(nrow(predict_coords)) - 1))
    valid_mask <- !is.na(idx)

    for (col in setdiff(selectors, "tmp_id")) {
      final_pred[idx[valid_mask], col] <- as.numeric(res_df[[col]][valid_mask])
    }
  } else {
    timestamp_message("Warning: Failed to process prediction coordinates (all dropped or GEE evaluation failed).")
  }

  # 5. Calculate Metrics
  if ("present" %in% names(final_pred)) {
    timestamp_message("  Calculating evaluation metrics...")
    
    # Check for Class Degeneracy (Must have both 0 and 1)
    u_vals <- unique(final_pred$present[!is.na(final_pred$present)])
    if (length(u_vals) < 2) {
        msg <- if (isTRUE(all(u_vals == 1))) "Presence-only" else "Absence-only"
        timestamp_message(sprintf("  ! Warning: Cannot calculate accuracy metrics (AUC/TSS). Evaluation set is %s.", msg))
        final_results$evaluation_error <- sprintf("Metric calculation skipped: validation set has %s samples.", msg)
    }
    
    metrics_list <- list()
    for (m in methods) {
      scores <- final_pred[[paste0("pred_", m)]]
      pos <- scores[final_pred$present == 1]
      neg <- scores[final_pred$present == 0]
      pos <- pos[!is.na(pos)]
      neg <- neg[!is.na(neg)]
      if (length(pos) > 0 && length(neg) > 0) metrics_list[[m]] <- calculate_classifier_metrics(pos, neg)
    }
    if (length(methods) > 1 && "pred_ensemble" %in% names(final_pred)) {
      scores <- final_pred$pred_ensemble
      pos <- scores[final_pred$present == 1]
      neg <- scores[final_pred$present == 0]
      pos <- pos[!is.na(pos)]
      neg <- neg[!is.na(neg)]
      if (length(pos) > 0 && length(neg) > 0) metrics_list[["ensemble"]] <- calculate_classifier_metrics(pos, neg)
    }
    final_results$metrics <- metrics_list

    # Live metric printing for ensemble
    if (!is.null(metrics_list[["ensemble"]])) {
      m <- metrics_list[["ensemble"]]
      timestamp_message(sprintf("  -> Evaluation (Ensemble): AUC=%.3f, TSS=%.3f", m$auc, m$tss))
    }
  }

  final_results$point_predictions <- final_pred
  return(final_results)
}
