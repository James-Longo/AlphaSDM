#' autoSDM: Automated Species Distribution Modeling
#'
#' This is the main entry point for the autoSDM pipeline.
#'
#' @param data A data frame formatted via `format_data()`. Must have standardized lowercase columns
#'   (longitude, latitude, year, present).
#' @param aoi Optional. Either a list with `lat`, `lon`, and `radius` (in meters), or a character string path to a polygon file.
#' @param output_dir Optional. Directory to save results. Defaults to the current working directory.
#' @param scale Optional. Resolution in meters for the final map. Defaults to 10.
#' @param python_path Optional. Path to Python executable. Auto-detected if not provided.
#' @param gee_project Optional. Google Cloud Project ID for Earth Engine.
#' @param cv Optional. Boolean whether to run spatial cross-validation. Defaults to FALSE.
#' @param predict_coords Optional. Data frame of coordinates to predict at.
#' @param methods Optional character vector of method names. Supported:
#'   Classifiers: `"rf"`, `"gbt"`, `"cart"`, `"maxent"`, `"svm"`.
#'   Reducers: `"centroid"`, `"ridge"`, `"linear"`, `"robust_linear"`.
#'   Defaults to all above methods.
#' @param ensemble Optional. Combine all methods into an ensemble map. Defaults to TRUE.
#' @param aoi_year Optional. Alpha Earth Mosaic year for map generation or background points. Defaults to NULL.
#' @param count Optional. Number of background points. Defaults to 10x presence points.
#' @param verbose_timer Optional. Boolean to display a real-time updating timer for each processing step. Defaults to FALSE.
#' @param n_cores Optional. Number of cores for parallel processing.
#' @param options Optional. A named list of advanced hyperparameters (e.g., `list(sampling_rate = 0.8, beta_multiplier = 2.0)`).
#'   Supported keys: `sampling_rate`, `loss_function`, `beta_multiplier`, `output_format`, `auto_feature`, `linear`, `quadratic`,
#'   `product`, `threshold`, `hinge`, `extrapolate`, `do_clamp`, `beta`, `seed`.
#' @return A list containing training data, models, and prediction results.
#' @export
autoSDM <- function(data,
                    aoi = NULL,
                    output_dir = "autoSDM_output",
                    scale = 10,
                    gee_project = NULL,
                    cv = 0,
                    predict_coords = NULL,
                    methods = NULL,
                    ensemble = TRUE,
                    aoi_year = 2017,
                    count = NULL,
                    n_trees = 100,
                    min_leaf_population = 5,
                    bag_fraction = 0.5,
                    shrinkage = 0.005,
                    max_nodes = NULL,
                    variables_per_split = NULL,
                    lambda_ = 0.1,
                    verbose_timer = FALSE,
                    options = list()) {

  # Set global timer option for nested calls
  original_timer <- getOption("autoSDM.verbose_timer")
  options(autoSDM.verbose_timer = verbose_timer)
  on.exit(options(autoSDM.verbose_timer = original_timer), add = TRUE)

  global_start <- log_time("Starting autoSDM pipeline")
  if (isTRUE(verbose_timer)) {
    message(sprintf("--- Training Configuration: n_trees=%d, shrinkage=%0.4f, min_leaf=%d ---",
                    n_trees, shrinkage, min_leaf_population))
  }
  # 1. Input Validation
  if (missing(data)) stop("Argument 'data' is required.")

  # Default AOI if not provided: Bounding box of data + buffer
  if (is.null(aoi) && is.null(predict_coords)) {
    message("No AOI provided. Calculating default bounding box from input data...")
    coords <- if (inherits(data, "sf")) sf::st_coordinates(data) else data[, c("longitude", "latitude")]

    min_lon <- min(coords[, 1], na.rm = TRUE)
    max_lon <- max(coords[, 1], na.rm = TRUE)
    min_lat <- min(coords[, 2], na.rm = TRUE)
    max_lat <- max(coords[, 2], na.rm = TRUE)

    lon_range <- max_lon - min_lon
    lat_range <- max_lat - min_lat

    # Center and radius approximation for CLI
    center_lon <- (min_lon + max_lon) / 2
    center_lat <- (min_lat + max_lat) / 2
    # Radius in meters approx (1 deg ~ 111km)
    # We use the larger dimension to cover the rectangle with a circle
    radius_m <- max(lon_range, lat_range) / 2 * 111000

    aoi <- list(lat = center_lat, lon = center_lon, radius = radius_m)
    message(sprintf("  Default AOI: %0.1f km radius around %0.4f, %0.4f", radius_m / 1000, center_lat, center_lon))
  }

  # Validate: need at least aoi or predict_coords
  if (is.null(aoi) && is.null(predict_coords)) {
    stop("You must provide either 'aoi' (for map generation) or 'predict_coords' (for point predictions), or both.")
  }
  # 1. Validate standardized column names
  required_cols <- c("longitude", "latitude", "year")
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    stop(sprintf(
      "Missing required columns: %s\nPlease use format_data() to standardize your data.",
      paste(missing, collapse = ", ")
    ))
  }


  # 2b. Check if data is Presence-Only (no absences)


  if (is.null(methods)) {
    methods <- c("centroid", "ridge", "linear", "robust_linear", "rf", "gbt", "cart", "maxent", "svm")
  }

  # Ensure "ensemble" is not accidentally listed in methods list
  methods <- setdiff(methods, "ensemble")

  if (length(methods) == 0) {
    stop("No valid matching methods found.")
  }

  # Build global tuning args
  tuning_args <- c(
    "--n-trees", as.character(as.integer(n_trees)),
    "--min-leaf-population", as.character(as.integer(min_leaf_population)),
    "--bag-fraction", as.character(bag_fraction),
    "--lambda", as.character(lambda_)
  )
  if (!is.null(shrinkage)) tuning_args <- c(tuning_args, "--shrinkage", as.character(shrinkage))
  if (!is.null(max_nodes)) tuning_args <- c(tuning_args, "--max-nodes", as.character(max_nodes))
  if (!is.null(variables_per_split)) tuning_args <- c(tuning_args, "--variables-per-split", as.character(variables_per_split))


  # Process advanced options
  if (length(options) > 0) {
    for (name in names(options)) {
      val <- options[[name]]
      cli_flag <- paste0("--", gsub("_", "-", name))
      tuning_args <- c(tuning_args, cli_flag, as.character(val))
    }
  }

  # If more than one method is specified, ensemble is typically expected, but we respect the explicit argument.

  # 3. GEE Configuration & Auth
  ensure_gee_authenticated(project = gee_project)
  ee <- reticulate::import("ee")

  # Resolve common parameters
  if (is.null(scale)) scale <- 10
  if (is.null(aoi_year)) aoi_year <- format(Sys.Date(), "%Y")

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # 5. Decide on Multi-Species
  is_multi_species <- "species" %in% names(data) && length(unique(data$species)) > 1
  working_data <- data

  if (!is_multi_species) {
    message("--- Single-Species Mode ---")

    # 6. GEE Geometry for AOI
    aoi_geom <- NULL
    if (is.list(aoi) && !is.null(aoi$lat)) {
      aoi_geom <- ee$Geometry$Point(c(as.numeric(aoi$lon), as.numeric(aoi$lat)))$buffer(as.numeric(aoi$radius))
    } else if (is.character(aoi) && file.exists(aoi)) {
      # Load from file using sf and convert to GEE
      aoi_sf <- sf::st_read(aoi, quiet = TRUE)
      aoi_geom <- rgee::sf_as_ee(aoi_sf)$geometry()
    } else if (!is.null(predict_coords)) {
      # Use bounding box of predict_coords
      bbox <- c(min(predict_coords$longitude), min(predict_coords$latitude),
                max(predict_coords$longitude), max(predict_coords$latitude))
      aoi_geom <- ee$Geometry$Rectangle(bbox)
    }

    # 7. Background Logic
    needs_bg <- any(methods %in% c("centroid", "ridge", "linear", "robust_linear", "rf", "gbt", "cart", "maxent", "svm")) &&
      !("present" %in% names(working_data) && any(working_data$present == 0))

    sampled_fc <- NULL
    if (needs_bg && !is.null(aoi_geom)) {
      message("--- Step: Generating Background Points (GEE-side) ---")
      n_bg <- if (!is.null(count)) count else (nrow(working_data) * 10)
      bg_fc <- get_background_embeddings_gee(aoi_geom, aoi_year, scale, n_bg)

      # Prepare presence points
      prep_pres <- prep_training_data_gee(working_data, class_property = "present", scale = scale)
      pres_fc <- prep_pres$fc

      sampled_fc <- pres_fc$merge(bg_fc)
    } else {
      prep <- prep_training_data_gee(working_data, class_property = "present", scale = scale)
      sampled_fc <- prep$fc
    }

    # 8. Training Loop
    message(sprintf("--- Training %d Method(s) via R-GEE ---", length(methods)))

    meta_results <- list()
    training_params <- list(
      numberOfTrees = n_trees,
      minLeafPopulation = min_leaf_population,
      bagFraction = bag_fraction,
      shrinkage = shrinkage,
      maxNodes = max_nodes,
      variablesPerSplit = variables_per_split,
      lambda_ = lambda_
    )

    # Add advanced options
    for (opt in names(options)) training_params[[opt]] <- options[[opt]]

    models <- list()
    for (m in methods) {
      msg <- sprintf("Training %s...", m)
      if (verbose_timer) message(msg)

      start_t <- Sys.time()
      models[[m]] <- train_gee_model(sampled_fc, m, params = training_params)
      end_t <- Sys.time()

      # Save metadata
      m_clean <- gsub("[^a-zA-Z0-9]", "", m)
      meta_path <- file.path(output_dir, paste0("model_", m_clean, ".json"))

      # For reducers, we store weights. For classifiers, we store the GEE ID if serialized?
      # Actually, let's keep the model objects in memory and just write a placeholder JSON.
      # (In the future, and for CRAN, we might want to serialize the GEE model ID)

      meta <- list(
        method = m,
        training_seconds = as.numeric(difftime(end_t, start_t, units = "secs")),
        is_classifier = models[[m]]$is_classifier
      )
      if (!is.null(models[[m]]$weights)) {
        meta$weights <- models[[m]]$weights
        meta$intercept <- models[[m]]$intercept
      }

      jsonlite::write_json(meta, meta_path, auto_unbox = TRUE)
      meta_results[[m]] <- meta_path
    }

    # 9. Map Generation (Individual)
    img_mosaic <- get_embedding_image(aoi_year, scale)

    if (!ensemble && !is.null(aoi_geom)) {
      for (m in methods) {
        message(sprintf("--- Generating Map: %s ---", m))
        pred_img <- predict_gee_map(models[[m]], img_mosaic)

        # Download or export
        # For now, let's assume local download if it's a small area
        m_clean <- gsub("[^a-zA-Z0-9]", "", m)
        tif_path <- file.path(output_dir, paste0(m_clean, ".tif"))

        # Use rgee::ee_as_raster for small, or task for large
        # We'll use a helper that decides based on size soon
        try(rgee::ee_as_raster(image = pred_img, region = aoi_geom, scale = scale, dsn = tif_path), silent = TRUE)
      }
    }

    # 10. Ensemble Logic
    final_results <- list(methods = methods, meta_files = meta_results)

    if (ensemble && !is.null(aoi_geom)) {
      message("--- Generating Ensemble Map ---")

      # consensus = mean of normalized scores
      # (This is simplified; a better ensemble would normalize first)
      pred_images <- lapply(models, function(mod) predict_gee_map(mod, img_mosaic))

      # Combine
      ensemble_img <- ee$ImageCollection(unname(pred_images))$mean()$rename("similarity")

      # Export
      tif_path <- file.path(output_dir, "ensemble.tif")
      try(rgee::ee_as_raster(image = ensemble_img, region = aoi_geom, scale = scale, dsn = tif_path), silent = TRUE)

      final_results$ensemble_map <- tif_path
    }

    # 11. Point Predictions
    if (!is.null(predict_coords)) {
      message("--- Predicting at specific coordinates (GEE-side) ---")
      final_results$point_predictions <- predict_at_coords_gee(
          coords_df = predict_coords,
          models = models,
          methods = methods,
          scale = scale,
          gee_project = gee_project
      )
    }

    log_time("autoSDM pipeline", global_start)
    return(final_results)
  } else {
    # MULTI-SPECIES WORKFLOW
    # -------------------------------------------------------------------------
    message(sprintf("--- Multi-Species Mode: Processing %d species ---", length(unique(working_data$species))))
    
    species_list <- unique(working_data$species)
    results_list <- list()
    
    for (sp in species_list) {
      message(sprintf("\n--- Processing Species: %s ---", sp))
      sp_data <- working_data[working_data$species == sp, ]
      sp_out <- file.path(output_dir, sp)
      
      # Handle predict_coords if scoped to species
      sp_predict <- predict_coords
      if (!is.null(predict_coords) && "species" %in% names(predict_coords)) {
        sp_predict <- predict_coords[predict_coords$species == sp, ]
      }
      
      # Recursive call for each species
      results_list[[sp]] <- autoSDM(
        data = sp_data,
        aoi = aoi,
        output_dir = sp_out,
        scale = scale,
        gee_project = gee_project,
        cv = cv,
        predict_coords = sp_predict,
        methods = methods,
        ensemble = ensemble,
        aoi_year = aoi_year,
        count = count,
        n_trees = n_trees,
        min_leaf_population = min_leaf_population,
        bag_fraction = bag_fraction,
        shrinkage = shrinkage,
        max_nodes = max_nodes,
        variables_per_split = variables_per_split,
        lambda_ = lambda_,
        verbose_timer = verbose_timer,
        options = options
      )
    }

    class(results_list) <- "autoSDM_multi"
    return(results_list)
  }
}
