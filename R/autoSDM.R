#' Generate SDM Map
#'
#' Trains models on provided data and generates spatial prediction maps (rasters) on GEE.
#'
#' @param data A data frame formatted via `format_data()`. Must have (longitude, latitude, year, present).
#'   If a 'species' column is present with multiple species, the function processes each separately.
#' @param aoi Either a list with `lat`, `lon`, and `radius` (in meters), or a path to a polygon file.
#' @param scale Resolution in meters for the final map. Defaults to 10.
#' @param output_dir Directory to save results. Defaults to the current working directory.
#' @param methods Character vector of method names. Supported: "rf", "gbt", "cart", "maxent", "svm", "ridge", "centroid".
#' @param ensemble Combine all methods into an ensemble map. Defaults to TRUE.
#' @param aoi_year Alpha Earth Mosaic year for map generation. Defaults to current year.
#' @param count Number of background points. Defaults to 10x presence points.
#' @param n_trees Number of trees for tree-based methods.
#' @param min_leaf_population Min leaf size for trees.
#' @param bag_fraction Bagging fraction.
#' @param shrinkage Learning rate for GBT.
#' @param lambda_ Penalty for Ridge.
#' @param polynomial Degree of polynomial terms for Ridge regression. 1 = linear, 2 = quadratic (default).
#' @param gee_project Optional. Google Cloud Project ID.
#' @param python_path Optional. Path to Python executable.
#' @param options Named list of advanced hyperparameters.
#' @export
generate_map <- function(data, aoi, scale = 10, output_dir = getwd(),
                        methods = NULL, ensemble = TRUE, aoi_year = NULL, count = NULL,
                        n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                        shrinkage = 0.005, max_nodes = NULL, variables_per_split = NULL, lambda_ = 0.0,
                        polynomial = 2L,
                        gee_project = NULL, python_path = NULL,
                        options = list()) {
  ensure_gee_authenticated(project = gee_project)
  ee <- reticulate::import("ee")

  # Multi-species check
  if ("species" %in% names(data) && length(unique(data$species)) > 1) {
    species_list <- unique(data$species)
    message(sprintf("--- Multi-Species Mode: Processing %d species ---", length(species_list)))
    results_list <- list()
    for (sp in species_list) {
      message(sprintf("\n--- Species: %s ---", sp))
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
  if (is.null(methods)) methods <- c("centroid", "ridge", "rf", "gbt", "cart", "maxent", "svm")
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

  train_res <- train_autoSDM_internal(
    data = data, methods = methods, aoi_geom = aoi_geom, scale = scale,
    aoi_year = aoi_year, count = count, training_params = training_params
  )

  models <- train_res$models
  final_results <- list(methods = methods, model_metadata = train_res$metadata)

  # 3. Map Generation
  img_mosaic <- get_embedding_image(aoi_year, scale)

  # Individual maps
  message("--- Generating Maps on GEE ---")
  for (m in methods) {
    pred_img <- predict_gee_map(models[[m]], img_mosaic)
    tif_path <- normalizePath(file.path(output_dir, paste0(gsub("[^a-zA-Z0-9]", "", m), ".tif")), mustWork = FALSE)
    message(sprintf("  Saving map to: %s", tif_path))
    try(rgee::ee_as_rast(image = pred_img, region = aoi_geom$bounds(), scale = scale, dsn = tif_path), silent = TRUE)
    final_results[[paste0(m, "_map")]] <- tif_path
  }

  # Ensemble map
  if (ensemble && length(methods) > 1) {
    message("  Calculating ensemble map...")
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
#' Trains models on provided data and evaluates them at specific coordinates to calculate accuracy metrics.
#'
#' @param data A data frame formatted via `format_data()`. Must have (longitude, latitude, year, present).
#' @param predict_coords Data frame of coordinates to evaluate at. Must have (longitude, latitude, year, present).
#' @param methods Character vector of method names.
#' @param scale Resolution in meters for sampling. Defaults to 10.
#' @param output_dir Directory to save results.
#' @param count Number of background points for training (if data is presence-only).
#' @param gee_project Optional. Google Cloud Project ID.
#' @param python_path Optional. Path to Python executable.
#' @param options Named list of advanced hyperparameters.
#' @export
evaluate_models <- function(data, predict_coords, methods = NULL, scale = 10, output_dir = getwd(),
                           aoi_year = NULL, count = NULL,
                           n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                           shrinkage = 0.005, max_nodes = NULL, variables_per_split = NULL, lambda_ = 0.0,
                           polynomial = 2L,
                           gee_project = NULL, python_path = NULL,
                           options = list()) {
  ensure_gee_authenticated(project = gee_project)
  ee <- reticulate::import("ee")

  # Multi-species check
  if ("species" %in% names(data) && length(unique(data$species)) > 1) {
    species_list <- unique(data$species)
    message(sprintf("--- Multi-Species Mode: Evaluating %d species ---", length(species_list)))
    results_list <- list()
    for (sp in species_list) {
      message(sprintf("\n--- Species: %s ---", sp))
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
  if (is.null(methods)) methods <- c("centroid", "ridge", "rf", "gbt", "cart", "maxent", "svm")
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

  train_res <- train_autoSDM_internal(
    data = data, methods = methods, aoi_geom = aoi_geom, scale = scale,
    aoi_year = aoi_year, count = count, training_params = training_params
  )

  models <- train_res$models
  final_results <- list(methods = methods, model_metadata = train_res$metadata)

  # 3. Coordinate Prediction
  message("--- Predicting at specific coordinates (Server-side) ---")
  upload_fc <- upload_points_to_gee(predict_coords)
  sampled_fc <- get_embeddings_at_fc(upload_fc, scale, properties = c("year", "tmp_id"))
  scored_fc <- predict_all_models_gee(sampled_fc, models)

  # Ensemble point scores
  if (length(methods) > 1) {
    message("  Calculating ensemble scores...")
    pred_cols <- paste0("pred_", methods)
    scored_fc <- scored_fc$map(function(f) {
      vals <- pred_cols %>% lapply(function(p) f$get(p))
      ensemble_score <- ee$List(vals)$reduce(ee$Reducer$mean())
      return(f$set("pred_ensemble", ensemble_score))
    })
  }

  # 4. Download Results
  selectors <- c("tmp_id", paste0("pred_", methods))
  if (length(methods) > 1) selectors <- c(selectors, "pred_ensemble")

  message(sprintf("  Downloading %d scores from GEE...", nrow(predict_coords)))
  result_table <- scored_fc$reduceColumns(
    reducer = ee$Reducer$toList()$`repeat`(length(selectors)),
    selectors = as.list(selectors)
  )$getInfo()

  res_df <- as.data.frame(do.call(cbind, result_table$list))
  colnames(res_df) <- selectors
  res_df$tmp_id <- as.numeric(res_df$tmp_id)
  res_df <- res_df[order(res_df$tmp_id), ]

  final_pred <- predict_coords
  for (col in setdiff(selectors, "tmp_id")) {
    final_pred[[col]] <- as.numeric(res_df[[col]])
  }

  # 5. Calculate Metrics
  if ("present" %in% names(final_pred)) {
    message("  Calculating evaluation metrics...")
    metrics_list <- list()
    for (m in methods) {
      scores <- final_pred[[paste0("pred_", m)]]
      pos <- scores[final_pred$present == 1]; neg <- scores[final_pred$present == 0]
      pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
      if (length(pos) > 0 && length(neg) > 0) metrics_list[[m]] <- calculate_classifier_metrics(pos, neg)
    }
    if (length(methods) > 1 && "pred_ensemble" %in% names(final_pred)) {
      scores <- final_pred$pred_ensemble
      pos <- scores[final_pred$present == 1]; neg <- scores[final_pred$present == 0]
      pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
      if (length(pos) > 0 && length(neg) > 0) metrics_list[["ensemble"]] <- calculate_classifier_metrics(pos, neg)
    }
    final_results$metrics <- metrics_list
  }

  final_results$point_predictions <- final_pred
  return(final_results)
}
