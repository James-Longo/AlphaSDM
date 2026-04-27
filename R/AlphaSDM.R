#' Internal Unified GEE Training Pipeline
#' @keywords internal
fit_gee_models <- function(train_df, methods, aoi_geom, scale, aoi_year, training_params) {
  ee <- reticulate::import("ee")

  # 1. Prepare Data and Sample Embeddings on GEE
  prep_res <- prep_training_data_gee(train_df, class_property = "present", scale = scale)
  sampled_fc <- prep_res$fc

  # 2. Train all requested models
  timestamp_message(sprintf("--- Training %d Method(s) via R-GEE ---", length(methods)))
  models <- list()
  for (m in methods) {
    timestamp_message(sprintf("  Training %s...", m))
    
    # Use method-specific params if available (from multi-tune), else global defaults
    m_params <- if (!is.null(training_params[[m]])) training_params[[m]] else training_params
    
    models[[m]] <- train_gee_model(sampled_fc, m, params = m_params)
  }

  return(list(
    models = models,
    metadata = list(
      n_presence = sum(train_df$present == 1),
      n_background = sum(train_df$present == 0),
      methods = methods,
      scale = scale
    )
  ))
}

#' Generate Background Points for AlphaSDM
#' @keywords internal
generate_background_points <- function(data, aoi_year, count = NULL, aoi = NULL) {
  ee <- reticulate::import("ee")
  
  # Coordinate Debug
  timestamp_message(sprintf("  -> Coord Range: Lon [%.2f, %.2f], Lat [%.2f, %.2f]", 
                            min(data$longitude), max(data$longitude), 
                            min(data$latitude), max(data$latitude)))
  
  if (!is.null(aoi)) {
    timestamp_message("--- Generating Background Points (Custom AOI) ---")
    bg_geoms <- list(aoi)
    bg_count_total <- if (is.null(count)) nrow(data) else count
    bg_counts <- list(bg_count_total)
  } else {
    timestamp_message("--- Generating Background Points (5-Cluster Spatial Strategy) ---")
    
    # 2a. Define Cluster regions
    if (nrow(data) >= 5) {
      bg_regions <- assign_spatial_folds(data, k = 5)$clusters
    } else {
      bg_regions <- list(st_as_sfc(st_bbox(st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326))))
    }
    
    # 2b. Convert SF geometries to GEE Geometries
    bg_geoms <- lapply(bg_regions, function(g) sf_as_ee(g))
    bg_count_total <- if (is.null(count)) nrow(data) else count
    bg_counts <- list(bg_count_total)
  }
  
  # 3. Quota-Safe Sampling (Oversample and filter for land points)
  timestamp_message(sprintf("  -> Target: %d land-based background points.", bg_count_total))
  
  # Combine all geometries into a single MultiPolygon for sampling
  # We use ee$FeatureCollection$randomPoints on a merged collection instead of MultiPolygon for safety
  bg_fcs_list <- lapply(bg_geoms, function(g) ee$FeatureCollection(ee$Feature(g)))
  combined_fc <- ee$FeatureCollection(bg_fcs_list)$flatten()
  
  get_valid_bg_df <- function(region_fc, total_to_request, target_count) {
    # We must add the 'year' property to the points so get_embeddings_at_fc knows which image to sample
    raw_fc <- ee$FeatureCollection$randomPoints(region_fc$geometry(), as.integer(total_to_request))
    
    # Add year property to every point
    raw_fc_with_year <- raw_fc$map(function(f) f$set("year", as.integer(aoi_year)))
    
    # Pass geometries = TRUE so we can retrieve the valid coordinates
    sampled_fc <- get_embeddings_at_fc(raw_fc_with_year, 100, geometries = TRUE)
    
    # Debug: Check sampled size
    actual_valid <- as.numeric(retry_curl_download(sampled_fc$size()$getInfo()))
    timestamp_message(sprintf("     [Debug] Requested %d, Found %d valid land points.", total_to_request, actual_valid))
    
    if (actual_valid == 0) return(NULL)
    
    info <- retry_curl_download(sampled_fc$limit(as.integer(target_count))$getInfo())
    
    do.call(rbind, lapply(info$features, function(f) {
      data.frame(
        longitude = f$geometry$coordinates[[1]],
        latitude = f$geometry$coordinates[[2]],
        year = aoi_year,
        present = 0
      )
    }))
  }
  
  # Initial attempt with 3x oversampling
  bg_df <- get_valid_bg_df(combined_fc, bg_count_total * 3, bg_count_total)
  
  # If we still don't have enough, retry with 10x
  if (is.null(bg_df) || nrow(bg_df) < bg_count_total) {
    timestamp_message("  -> Insufficient land points in first pass, retrying with 10x oversampling...")
    bg_df <- get_valid_bg_df(combined_fc, bg_count_total * 10, bg_count_total)
  }
  
  if (is.null(bg_df)) stop("Could not find any valid land points for background sampling.")
  
  timestamp_message(sprintf("  -> Successfully secured %d valid background points.", nrow(bg_df)))
  return(rbind(data[, c("longitude", "latitude", "year", "present")], bg_df))
}

#' Generate SDM Map
#' @export
generate_map <- function(data, aoi, scale = 10, output_dir = getwd(),
                         methods = NULL, ensemble = TRUE, aoi_year = NULL, count = NULL,
                         n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                         shrinkage = 0.005, max_nodes = 6L, variables_per_split = NULL,
                         gee_project = NULL, python_path = NULL,
                         options = list()) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)
  ee <- reticulate::import("ee")

  if (is.null(aoi_year)) aoi_year <- 2023
  if (is.null(methods)) methods <- c("similarity", "rf", "gbt", "maxent", "svm")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. Prepare AOI Geometry
  aoi_geom <- NULL
  if (is.list(aoi) && !is.null(aoi$lat)) {
    aoi_geom <- ee$Geometry$Point(c(as.numeric(aoi$lon), as.numeric(aoi$lat)))$buffer(as.numeric(aoi$radius))
  } else if (is.character(aoi) && file.exists(aoi)) {
    aoi_sf <- sf::st_read(aoi, quiet = TRUE)
    aoi_geom <- rgee::sf_as_ee(aoi_sf)$geometry()
  }

  # 2. Prepare Training Data with Background
  if (all(data$present == 1)) {
    data <- generate_background_points(data, aoi_year, count, aoi = aoi_geom)
  }

  # 3. Method-Specific Optimized Defaults (Applied if not explicitly tuned)
  training_params <- list(
    numberOfTrees = n_trees, 
    minLeafPopulation = min_leaf_population,
    bagFraction = bag_fraction, 
    shrinkage = shrinkage,
    maxNodes = max_nodes, 
    variablesPerSplit = variables_per_split
  )

  # Override GBT with 150 trees if using global default of 100
  if ("gbt" %in% methods && n_trees == 100L) {
    training_params$gbt <- training_params
    training_params$gbt$numberOfTrees <- 150L
  }
  
  # Override RF with 250 trees if using global default of 100
  if ("rf" %in% methods && n_trees == 100L) {
    training_params$rf <- training_params
    training_params$rf$numberOfTrees <- 250L
  }

  # 4. Unified Training
  train_res <- fit_gee_models(data, methods, aoi_geom, scale, aoi_year, training_params)

  # 5. Map Generation
  img_mosaic <- get_embedding_image(aoi_year, scale)
  final_results <- list(methods = methods, model_metadata = train_res$metadata)
  for (m in methods) {
    pred_img <- predict_gee_map(train_res$models[[m]], img_mosaic)
    tif_path <- file.path(output_dir, paste0(m, ".tif"))
    try(rgee::ee_as_rast(image = pred_img, region = aoi_geom$bounds(), scale = scale, dsn = tif_path), silent = TRUE)
    final_results[[paste0(m, "_map")]] <- tif_path
  }
  return(final_results)
}

#' Evaluate SDM Models
#' @export
evaluate_models <- function(data, predict_coords, scale = 10, output_dir = getwd(),
                            methods = NULL, aoi_year = NULL, count = NULL,
                            n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                            shrinkage = 0.005, max_nodes = 6L, variables_per_split = NULL,
                            gee_project = NULL, python_path = NULL,
                            options = list()) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)
  ee <- reticulate::import("ee")

  if (is.null(aoi_year)) aoi_year <- 2023
  if (is.null(methods)) methods <- c("similarity", "rf", "gbt", "maxent", "svm")

  # 1. Prepare AOI Geometry (Bounding Box of predict_coords)
  bbox <- c(min(predict_coords$longitude), min(predict_coords$latitude), max(predict_coords$longitude), max(predict_coords$latitude))
  aoi_geom <- ee$Geometry$Rectangle(bbox)

  # 2. Background Points (Honest: Constrained to Training Data Extent)
  if (all(data$present == 1)) {
    data <- generate_background_points(data, aoi_year, count)
  }

  # 3. Method-Specific Optimized Defaults
  training_params <- list(
    numberOfTrees = n_trees, 
    minLeafPopulation = min_leaf_population,
    bagFraction = bag_fraction, 
    shrinkage = shrinkage,
    maxNodes = max_nodes, 
    variablesPerSplit = variables_per_split
  )

  # Override GBT with 150 trees if using global default of 100
  if ("gbt" %in% methods && n_trees == 100L) {
    training_params$gbt <- training_params
    training_params$gbt$numberOfTrees <- 150L
  }
  
  # Override RF with 250 trees if using global default of 100
  if ("rf" %in% methods && n_trees == 100L) {
    training_params$rf <- training_params
    training_params$rf$numberOfTrees <- 250L
  }

  # 4. Unified Training
  train_res <- fit_gee_models(data, methods, aoi_geom, scale, aoi_year, training_params)

  # 3. Prediction Pipeline
  if (!"year" %in% names(predict_coords)) {
    predict_coords$year <- aoi_year
  }

  timestamp_message(sprintf("--- Predicting at %d coordinates (Server-side, Batched) ---", nrow(predict_coords)))
  chunk_size <- if (!is.null(options$batch_size)) as.integer(options$batch_size) else 5000L
  batch_starts <- seq(1, nrow(predict_coords), by = chunk_size)

  process_batch <- function(idx) {
    try({
      i <- batch_starts[idx]
      end_idx <- min(i + chunk_size - 1, nrow(predict_coords))
      chunk_df <- predict_coords[i:end_idx, , drop = FALSE]
      upload_fc <- upload_points_to_gee(chunk_df)
      sampled_fc <- get_embeddings_at_fc(upload_fc, scale, properties = c("year"))
  
      # Run all model predictions (Classifiers & Reducers)
      sampled_fc <- predict_all_models_gee(sampled_fc, train_res$models)
      return(retry_curl_download(sampled_fc$getInfo()))
    }, silent = TRUE)
  }

  # Use future_lapply for parallelism (user-controlled)
  all_res_list <- future_lapply(seq_along(batch_starts), process_batch, future.seed = TRUE)

  # 4. Metrics & Result Collection
  final_pred <- predict_coords
  for (m in methods) {
    col <- paste0("pred_", m)
    
    # Safely extract predictions, skipping failed batches
    all_preds <- unlist(lapply(all_res_list, function(res) {
      if (inherits(res, "try-error") || is.null(res) || !"features" %in% names(res)) {
        return(rep(NA, chunk_size)) # Placeholder for failed batch (approximate)
      }
      sapply(res$features, function(f) {
        val <- f$properties[[col]]
        if (is.null(val)) NA else as.numeric(val)
      })
    }))
    
    # Ensure length matches (pad if needed)
    if (length(all_preds) < nrow(final_pred)) {
      all_preds <- c(all_preds, rep(NA, nrow(final_pred) - length(all_preds)))
    } else if (length(all_preds) > nrow(final_pred)) {
      all_preds <- all_preds[1:nrow(final_pred)]
    }
    
    final_pred[[col]] <- all_preds
  }

  final_results <- list(methods = methods, metrics = list(), point_predictions = final_pred)
  if ("present" %in% names(final_pred)) {
    for (m in methods) {
      scores <- final_pred[[paste0("pred_", m)]]
      pos <- scores[final_pred$present == 1]
      neg <- scores[final_pred$present == 0]
      final_results$metrics[[m]] <- calculate_classifier_metrics(pos, neg)
    }
  }
  return(final_results)
}
