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
autoSDM <- function(data, aoi = NULL, output_dir = getwd(), scale = NULL, python_path = NULL,
                    gee_project = NULL, cv = FALSE, predict_coords = NULL,
                    methods = NULL, ensemble = TRUE, aoi_year = NULL, count = NULL,
                    n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                    shrinkage = 0.005, max_nodes = NULL, variables_per_split = NULL, lambda_ = 0.1,
                    n_cores = NULL, verbose_timer = FALSE, options = list()) {
  # 1. Input Validation
  if (missing(data)) stop("Argument 'data' is required.")
  if (verbose_timer) message("--- autoSDM: Real-time processing feedback enabled ---")

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

  message(sprintf("--- Training Configuration: n_trees=%d, shrinkage=%.4f, min_leaf=%d ---", n_trees, shrinkage, min_leaf_population))

  # Process advanced options
  if (length(options) > 0) {
    for (name in names(options)) {
      val <- options[[name]]
      cli_flag <- paste0("--", gsub("_", "-", name))
      tuning_args <- c(tuning_args, cli_flag, as.character(val))
    }
  }

  # If more than one method is specified, ensemble is typically expected, but we respect the explicit argument.

  # 3. Python Configuration
  # Check for virtualenv and initialize dependencies
  python_path_detected <- ensure_autoSDM_dependencies()
  python_path <- if (!is.null(python_path)) python_path else python_path_detected
  python_path <- resolve_python_path(python_path)

  if (is.null(python_path) || !file.exists(python_path)) {
    stop("Could not find a valid Python environment.\nPlease ensure Python is installed and detected by `reticulate::py_config()`, or provide the `python_path` argument explicitly.")
  }

  # Locate the python source directory (inst/python)
  pkg_py_path <- ""
  if (file.exists(file.path(getwd(), "inst", "python"))) {
    pkg_py_path <- file.path(getwd(), "inst", "python")
  } else {
    pkg_py_path <- system.file("python", package = "autoSDM")
  }

  if (pkg_py_path != "") {
    Sys.setenv(PYTHONPATH = pkg_py_path)
    message(sprintf("Added to PYTHONPATH: %s", pkg_py_path))
  }

  # 4. Check GEE Readiness
  message(sprintf("Using Python: %s", python_path))
  ensure_gee_authenticated(project = gee_project)

  # If project was not provided, try to load it from the config that ensure_gee_authenticated might have just used/created
  if (is.null(gee_project)) {
    config_file <- file.path(Sys.getenv("HOME"), ".config", "autoSDM", "config.json")
    if (file.exists(config_file)) {
      try(
        {
          conf <- jsonlite::fromJSON(config_file)
          if (!is.null(conf$gee_project)) {
            gee_project <- conf$gee_project
            message("Loaded GEE Project from config: ", gee_project)
          }
        },
        silent = TRUE
      )
    }
  }

  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Temporary directory for inter-process CV cache (never written to output_dir)
  cv_cache_dir <- tempfile(pattern = "autoSDM_cv_cache_")
  dir.create(cv_cache_dir, showWarnings = FALSE)
  cv_cache_arg <- c("--cv-cache-dir", shQuote(cv_cache_dir))

  # 5. Decide on Multi-Species
  is_multi_species <- "species" %in% names(data) && length(unique(data$species)) > 1

  # Use the raw input data as our base; analyze/predict CLI will sample GEE as needed.
  working_data <- data

  # 5. Prepare AOI arguments
  aoi_arg <- NULL
  if (is.list(aoi) && !is.null(aoi$lat)) {
    aoi_arg <- c("--lat", aoi$lat, "--lon", aoi$lon, "--radius", aoi$radius)
  } else if (is.character(aoi)) {
    aoi_arg <- c("--aoi-path", shQuote(aoi))
  }

  # 6. Extraction/Preparation
  # (Standardizing CSV output for Python)
  extract_csv <- tempfile(fileext = ".csv")
  # working_data is our raw input (coords only)
  if (requireNamespace("vroom", quietly = TRUE)) {
    vroom::vroom_write(working_data, extract_csv, delim = ",")
  } else {
    write.csv(working_data, extract_csv, row.names = FALSE)
  }

  if (!is_multi_species) {
    message("--- Single-Species Mode ---")
    proj_arg <- if (!is.null(gee_project)) c("--project", shQuote(gee_project)) else NULL

    meta_files <- list()

    # Shared Background logic for Presence-Only
    # If any model needs background (centroid/ridge) and we don't have absences yet.
    needs_bg <- any(methods %in% c("centroid", "ridge")) &&
      !("present" %in% names(working_data) && any(working_data$present == 0))

    if (needs_bg && (!is.null(aoi) || !is.null(aoi_arg))) {
      message("--- Step: Generating Shared Background Points (10:1) ---")
      bg_csv <- tempfile(fileext = ".csv")
      n_bg <- nrow(working_data) * 10

      .run_command_with_timer(python_path, args = c(
        "-m", "autoSDM.cli", "background",
        "--output", shQuote(bg_csv),
        "--count", if (!is.null(count)) count else n_bg,
        "--scale", scale,
        if (!is.null(aoi_year)) c("--year", aoi_year) else NULL,
        proj_arg, aoi_arg
      ), label = "Generating background points", active = verbose_timer)

      if (file.exists(bg_csv)) {
        bg_data <- read.csv(bg_csv)
        if (!"present" %in% names(bg_data)) bg_data$present <- 0
        if (!"present" %in% names(working_data)) working_data$present <- 1

        # Ensure column alignment
        common_cols <- intersect(names(working_data), names(bg_data))
        working_data <- rbind(working_data[, common_cols], bg_data[, common_cols])

        # Update standardized CSV with background points
        if (requireNamespace("vroom", quietly = TRUE)) {
          vroom::vroom_write(working_data, extract_csv, delim = ",")
        } else {
          write.csv(working_data, extract_csv, row.names = FALSE)
        }
      }
    }

    # Track execution time
    execution_times <- list()

    message(sprintf("--- Dispatching %d Method(s) via Ensemble (GEE-side) ---", length(methods)))

    # We call the CLI ONCE with all methods to "keep all the data in GEE" 
    # and avoid the 10MB payload limit associated with uploading embeddings.
    methods_csv <- paste(methods, collapse = ",")
    
    start_time_all <- Sys.time()
    .run_command_with_timer(python_path, args = c(
      "-m", "autoSDM.cli", "analyze",
      "--input", shQuote(extract_csv),
      "--output", shQuote(file.path(output_dir, "model.csv")), # Will generate model_method.csv for each
      "--method", methods_csv,
      "--scale", scale,
      if (!is.null(aoi_year)) c("--year", aoi_year) else NULL,
      if (!is.null(count)) c("--count", count) else NULL,
      tuning_args, cv_cache_arg, proj_arg
    ), label = "Training GEE models (Shared Context)", active = verbose_timer)
    end_time_all <- Sys.time()

    # Reconstruct parallel_results for compatibility with the rest of the R script
    parallel_results <- lapply(methods, function(m) {
      m_clean <- gsub("[^a-zA-Z0-9]", "", m)
      # CLI now outputs model_method.csv and model_method.json
      out_json <- file.path(output_dir, paste0("model_", m_clean, ".json"))
      
      ret <- list(method = m, meta_file = NULL, training_seconds = as.numeric(difftime(end_time_all, start_time_all, units = "secs")) / length(methods))
      if (file.exists(out_json)) {
        ret$meta_file <- out_json
      }
      return(ret)
    })

    # Reconstruct meta_files and execution_times
    for (res in parallel_results) {
      m <- res$method
      m_clean <- gsub("[^a-zA-Z0-9]", "", m)
      out_json <- res$meta_file

      if (!is.null(out_json)) {
        meta_files[[m]] <- out_json
        
        # 7. Generate Individual Map for this method (if requested and not in ensemble mode)
        if (!is.null(aoi) && !ensemble) {
           message(sprintf("--- Step: Generating Individual Map (%s) ---", m))
           .run_command_with_timer(python_path, args = c(
             "-m", "autoSDM.cli", "extrapolate",
             "--meta", shQuote(out_json),
             "--scale", scale,
             "--prefix", m_clean,
             if (!is.null(aoi_year)) c("--year", aoi_year) else NULL,
             proj_arg, aoi_arg
           ), label = paste("Mapping", m_clean, "predictions"), active = verbose_timer)
        }
      }
      
      if (is.null(execution_times[[m]])) execution_times[[m]] <- list()
      execution_times[[m]]$training_seconds <- res$training_seconds
    }

    # 8. Predict at specific coordinates (if provided)
    # 8. Predict at specific coordinates (if provided)
    if (!is.null(predict_coords)) {
      message("--- Step: Predicting at specific coordinates (GEE-side orchestration) ---")
      
      point_preds <- predict_coords
      submodel_test_metrics <- list()

      # Collect all valid meta files
      valid_metas <- as.character(unlist(meta_files))
      valid_metas <- valid_metas[file.exists(valid_metas)]

      if (length(valid_metas) > 0) {
        pred_start_time <- Sys.time()
        # Call predict_at_coords once for all methods
        preds_all <- predict_at_coords(
          predict_coords, 
          analysis_meta_paths = valid_metas, 
          scale = scale, 
          aoi_year = aoi_year, 
          python_path = python_path, 
          gee_project = gee_project, 
          verbose_timer = verbose_timer
        )
        pred_end_time <- Sys.time()
        
        point_preds <- preds_all$data
        point_preds$similarity <- 0
        method_count <- 0
        submodel_test_metrics <- list()

        # Update point_preds with normalized similarities and individual results
        for (m in names(meta_files)) {
          m_clean <- gsub("[^a-zA-Z0-9]", "", m)
          col_name <- if (length(valid_metas) > 1) paste0("similarity_", m) else "similarity"
          
          if (!col_name %in% names(point_preds)) next

          # Load meta to get similarity range for normalization
          m_meta <- jsonlite::fromJSON(meta_files[[m]])
          s_range <- m_meta$similarity_range
          
          raw_sim <- point_preds[[col_name]]
          if (!is.null(s_range) && (s_range[2] - s_range[1] > 1e-9)) {
            norm_sim <- (raw_sim - s_range[1]) / (s_range[2] - s_range[1])
          } else {
            norm_sim <- raw_sim
          }
          
          # Store normalized back in the df using method name as column
          point_preds[[m]] <- norm_sim
          point_preds$similarity <- point_preds$similarity + norm_sim
          method_count <- method_count + 1
          
          # Handle metrics
          if (!is.null(preds_all$metrics[[m]])) {
            submodel_test_metrics[[m]] <- preds_all$metrics[[m]]$testing
          }
          
          if (is.null(execution_times[[m]])) execution_times[[m]] <- list()
          execution_times[[m]]$prediction_seconds <- as.numeric(difftime(pred_end_time, pred_start_time, units = "secs")) / length(valid_metas)
        }

        if (method_count > 0) {
          point_preds$similarity <- point_preds$similarity / method_count
        }

        # Final normalization of the ensemble product to 0-1
        e_min <- min(point_preds$similarity, na.rm = TRUE)
        e_max <- max(point_preds$similarity, na.rm = TRUE)
        if (e_max - e_min > 1e-9) {
          point_preds$similarity <- (point_preds$similarity - e_min) / (e_max - e_min)
        }
      }
    }

    # 9. Generate Ensemble Map (if requested)
    if (!is.null(aoi) && ensemble) {
      message("--- Step: Generating Ensemble Map (consensus of all methods) ---")

      ensemble_results_json <- file.path(output_dir, "ensemble_results.json")

      # Build args with all meta files
      ensemble_args <- c("-m", "autoSDM.cli", "ensemble", "--input", shQuote(extract_csv), "--output", shQuote(ensemble_results_json))
      for (mf in meta_files) {
        if (file.exists(mf)) ensemble_args <- c(ensemble_args, "--meta", shQuote(mf))
      }
      ensemble_args <- c(ensemble_args, "--scale", scale, "--prefix", "ensemble")
      if (!is.null(aoi_year)) ensemble_args <- c(ensemble_args, "--year", aoi_year)

      if (!is.null(gee_project)) ensemble_args <- c(ensemble_args, "--project", shQuote(gee_project))

      if (is.list(aoi) && !is.null(aoi$lat)) {
        ensemble_args <- c(ensemble_args, "--lat", aoi$lat, "--lon", aoi$lon, "--radius", aoi$radius)
      } else if (is.character(aoi)) {
        ensemble_args <- c(ensemble_args, "--aoi-path", shQuote(aoi))
      }

      # Run CV on ensemble formulation
      cv_arg_ensemble <- if (cv) {
        c("--cv", "--train-methods", paste(methods, collapse = ","), "--eval-methods", "ensemble")
      } else {
        NULL
      }

      if (!is.null(cv_arg_ensemble)) ensemble_args <- c(ensemble_args, cv_arg_ensemble, cv_cache_arg)

      status <- .run_command_with_timer(python_path, args = ensemble_args, label = "Creating ensemble map", active = verbose_timer)
      if (status != 0) stop("Ensemble extrapolation failed.")
      final_results <- jsonlite::fromJSON(ensemble_results_json)
    } else {
      # For single results, we return the last json generated or empty list
      # The individual maps are already handled in the loop.
      final_results <- list()
      if (length(meta_files) > 0) {
        final_results <- jsonlite::fromJSON(meta_files[[1]])
      }
    }

    if (!is.null(predict_coords)) {
      final_results$point_predictions <- point_preds
      
      # If we have an ensemble result, calculate its testing metrics if labels are present
      if (ensemble && !is.null(point_preds$present) && any(!is.na(point_preds$present))) {
          # Use the testing metrics from submodels if it's a single method, 
          # but for ensemble we might need a separate calculation.
          # For now, if it's a single method, we just take its test_metrics.
          if (length(methods) == 1 && !is.null(submodel_test_metrics[[methods[1]]])) {
              final_results$metrics$testing <- submodel_test_metrics[[methods[1]]]
          } else if (length(methods) > 1) {
              # OPTIONAL: Add R-side ensemble AUC/CBI calculation here.
              # For now, we'll mark it as available per sub-model in execution_times or similar,
              # but let's try to at least provide the first one if it's representative or empty.
          }
      } else if (!ensemble && length(methods) == 1) {
          final_results$metrics$testing <- submodel_test_metrics[[methods[1]]]
      }
    }
    
    if (length(execution_times) > 0) final_results$execution_times <- execution_times

    message("autoSDM pipeline complete!")
    return(final_results)
  } else {
    # MULTI-SPECIES WORKFLOW
    # -------------------------------------------------------------------------
    message(sprintf("--- Multi-Species Mode: Processing %d species ---", length(unique(working_data$species))))
    message("Orchestrating GEE server-side parallelization...")


    cv_arg <- if (cv) "--cv" else NULL
    proj_arg <- if (!is.null(gee_project)) c("--project", shQuote(gee_project)) else NULL

    aoi_arg <- NULL
    if (is.list(aoi) && !is.null(aoi$lat)) {
      aoi_arg <- c("--lat", aoi$lat, "--lon", aoi$lon, "--radius", aoi$radius)
    } else if (is.character(aoi)) {
      aoi_arg <- c("--aoi-path", shQuote(aoi))
    }

    # Internal multi-species ensemble pipeline in Python
    # This invokes a single Python process that loops over species and sends requests to GEE.
    multi_args <- c(
      "-m", "autoSDM.cli", "analyze",
      "--input", shQuote(extract_csv),
      "--output", shQuote(output_dir),
      "--method", "ensemble",
      "--species-col", "species",
      "--scale", scale,
      if (!is.null(aoi_year)) c("--year", aoi_year) else NULL,
      cv_arg,
      proj_arg,
      aoi_arg,
      tuning_args
    )

    status <- .run_command_with_timer(python_path, args = multi_args, label = "Running multi-species batch processing", active = verbose_timer)
    if (status != 0) stop("Multi-species analysis failed.")

    # Load and combine results from species-specific directories
    species_list <- unique(working_data$species)
    results_list <- list()
    for (sp in species_list) {
      res_path <- file.path(output_dir, sp, "results.json")
      if (file.exists(res_path)) {
        results_list[[sp]] <- jsonlite::fromJSON(res_path)
      }
    }

    class(results_list) <- "autoSDM_multi"
  }
}
