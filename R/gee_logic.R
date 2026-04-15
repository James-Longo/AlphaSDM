#' Internal helper for retrying GEE operations with exponential backoff
#'
#' @keywords internal
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
      timestamp_message(sprintf("  GEE rate limit/timeout hit. Retrying in %.2f seconds (Attempt %d/%d)...", delay, i, max_retries))
      Sys.sleep(delay)
    } else {
      stop(res)
    }
  }
}

get_embedding_image <- function(year, scale = 10) {
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
#' @keywords internal
get_embeddings_at_fc_raw <- function(fc, scale, properties = NULL, geometries = FALSE, years = NULL) {
  ee <- reticulate::import("ee")

  if (is.null(years)) {
    years <- fc$aggregate_array("year")$distinct()$getInfo()
  }

  sampled_fcs <- list()
  for (yr in years) {
    yr_fc <- fc$filter(ee$Filter$eq("year", as.integer(yr)))
    yr_img <- get_embedding_image(yr, scale)

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

#' Prepare Training Data for GEE
#' Prepare Training Data on GEE (Pure Server-Side)
#'
#' @param df Data frame with longitude, latitude, year, and class_property
#' @param class_property Column name for presence/absence
#' @param scale Resolution in meters
#' @keywords internal
prep_training_data_gee <- function(df, class_property = "present", scale = 10) {
  ee <- reticulate::import("ee")

  # 1. Clean data
  cols_to_keep <- c("longitude", "latitude", "year")
  if (!is.null(class_property)) cols_to_keep <- c(cols_to_keep, class_property)

  df_clean <- df[complete.cases(df[, cols_to_keep]), ]
  if (nrow(df_clean) == 0) stop("No valid training data remaining after dropping NAs.")

  # 2. Upload to GEE
  # For training, we group by class/year for efficiency
  timestamp_message("\n--- Uploading User Data ---")
  timestamp_message(sprintf("  -> Transferring %d provided coordinates to GEE...", nrow(df_clean)))

  group_cols <- "year"
  if (!is.null(class_property)) group_cols <- c(group_cols, class_property)

  group_fac <- if (length(group_cols) > 1) {
    interaction(df_clean[, group_cols], drop = TRUE)
  } else {
    df_clean[[group_cols]]
  }
  groups <- split(df_clean, group_fac)

  gee_features <- list()
  for (grp_name in names(groups)) {
    grp <- groups[[grp_name]]
    if (nrow(grp) == 0) next

    yr <- grp$year[1]
    cls_val <- if (!is.null(class_property)) grp[[class_property]][1] else NULL

    # Props for this group
    props <- list(year = as.integer(yr))
    if (!is.null(class_property)) props[[class_property]] <- cls_val

    # Convert points to MultiPoint
    coords_mat <- as.matrix(grp[, c("longitude", "latitude")])
    chunk_size <- 1000
    for (i in seq(1, nrow(coords_mat), by = chunk_size)) {
      end <- min(i + chunk_size - 1, nrow(coords_mat))
      sub_coords <- coords_mat[i:end, , drop = FALSE]
      coord_list <- lapply(seq_len(nrow(sub_coords)), function(j) as.list(unname(sub_coords[j, ])))
      geom <- ee$Geometry$MultiPoint(coord_list)
      gee_features <- c(gee_features, list(ee$Feature(geom, props)))
    }
  }

  compressed_fc <- ee$FeatureCollection(gee_features)

  # 2.5 Unpack MultiPoints into individual Point features on GEE
  # This ensures sampleRegions treats every coordinate as a unique point.
  reticulate::py_run_string("
import ee
def unpack_simple_multipoint(f):
    f = ee.Feature(f)
    coords = f.geometry().coordinates()
    yr = f.get('year')
    cls = f.get('present')
    return coords.map(
        lambda pt: ee.Feature(ee.Geometry.Point(pt), {'year': yr, 'present': cls})
    )
  ")

  unpack_func <- reticulate::py$unpack_simple_multipoint
  num_feats <- as.integer(length(gee_features))
  upload_fc <- compressed_fc$toList(num_feats)$map(unpack_func)$flatten()
  upload_fc <- ee$FeatureCollection(upload_fc)

  # Logging size
  timestamp_message(sprintf("GEE FeatureCollection created & unpacked: %d Points.", as.integer(upload_fc$size()$getInfo())))

  # 3. Sample embeddings
  sample_props <- group_cols
  sampled_fc <- get_embeddings_at_fc(upload_fc, scale, properties = sample_props)

  # Log valid points remaining
  final_count <- as.numeric(sampled_fc$size()$getInfo())
  dropped <- nrow(df_clean) - final_count
  if (dropped > 0) {
    timestamp_message(sprintf("  -> Discarded %d presence points due to missing embeddings (e.g., over water).", dropped))

    # Attempt to retrieve the discarded coordinates for diagnostic investigation
    tryCatch(
      {
        # Use an inverted join to find features that were in the upload but not in the result
        # This works by finding features in 'upload_fc' whose geometry doesn't match anything in 'sampled_fc'
        dist_filter <- ee$Filter$withinDistance(
          distance = scale, # Small spatial tolerance
          leftField = ".geo",
          rightField = ".geo"
        )

        inverted_join <- ee$Join$inverted()
        orphans_fc <- inverted_join$apply(upload_fc, sampled_fc, dist_filter)

        # Limit to first 500 orphans to avoid GEE memory overflow during download
        orphans_info <- orphans_fc$limit(500L)$getInfo()$features

        if (length(orphans_info) > 0) {
          orphan_rows <- list()
          for (f in orphans_info) {
            # MultiPoint coordinates is a list of [lon, lat] pairs
            coords_list <- f$geometry$coordinates
            for (pt in coords_list) {
              orphan_rows <- c(orphan_rows, list(data.frame(
                longitude = pt[[1]],
                latitude = pt[[2]],
                year = f$properties$year,
                species = if ("species" %in% names(df_clean)) df_clean$species[1] else "unknown"
              )))
            }
          }
          orphan_df <- do.call(rbind, orphan_rows)

          # Log to shared CSV
          log_path <- "discarded_embedding_points.csv"
          write.table(orphan_df,
            file = log_path, sep = ",", append = file.exists(log_path),
            col.names = !file.exists(log_path), row.names = FALSE
          )
          timestamp_message(sprintf("  -> Logged %d orphaned coordinates to %s", nrow(orphan_df), log_path))
        }
      },
      error = function(e) {
        timestamp_message(sprintf("  ! Failed to log orphaned coordinates: %s", e$message))
      }
    )
  }
  timestamp_message(sprintf("  -> Final valid presence points: %d", final_count))

  return(list(
    fc = sampled_fc,
    df_clean = df_clean,
    class_property = class_property
  ))
}

#' Upload Individual Points to GEE with IDs
#'
#' @param df Data frame with longitude, latitude, year
#' @param chunk_size Points per upload chunk
#' @param id_offset Numeric offset for tmp_id to support outer batching
#' @return GEE FeatureCollection
#' @keywords internal
upload_points_to_gee <- function(df, chunk_size = 2000, id_offset = 0) {
  ee <- reticulate::import("ee")
  timestamp_message(sprintf("Uploading %d coordinates to GEE via efficient MultiPoint compression...", nrow(df)))

  df$tmp_id <- as.integer(id_offset + seq_len(nrow(df)) - 1)
  groups <- split(df, df$year)

  all_features <- list()
  for (yr_str in names(groups)) {
    grp <- groups[[yr_str]]
    yr <- as.integer(yr_str)

    for (i in seq(1, nrow(grp), by = chunk_size)) {
      end <- min(i + chunk_size - 1, nrow(grp))
      chunk_df <- grp[i:end, , drop = FALSE]

      coords_list <- lapply(seq_len(nrow(chunk_df)), function(j) {
        list(as.numeric(chunk_df$longitude[j]), as.numeric(chunk_df$latitude[j]))
      })
      ids_list <- as.list(as.integer(chunk_df$tmp_id))

      geom <- ee$Geometry$MultiPoint(coords_list)
      props <- list(
        year = yr,
        tmp_ids = ee$List(ids_list),
        lons = ee$List(as.list(as.numeric(chunk_df$longitude))),
        lats = ee$List(as.list(as.numeric(chunk_df$latitude)))
      )

      all_features <- c(all_features, list(ee$Feature(geom, props)))
    }
  }

  compressed_fc <- ee$FeatureCollection(unname(all_features))
  num_feats <- as.integer(length(all_features))

  # Unpack cleanly on remote server using pure Python definition to bypass R scoping limits
  reticulate::py_run_string("
import ee
def unpack_compressed_geom(f):
    f = ee.Feature(f)
    coords = f.geometry().coordinates()
    ids = ee.List(f.get('tmp_ids'))
    lons = ee.List(f.get('lons'))
    lats = ee.List(f.get('lats'))
    yr = f.get('year')
    return coords.zip(ids).zip(lons).zip(lats).map(
        lambda triple: ee.Feature(
            ee.Geometry.Point(ee.List(ee.List(ee.List(triple).get(0)).get(0)).get(0)),
            {
                'year': yr,
                'tmp_id': ee.List(ee.List(ee.List(triple).get(0)).get(0)).get(1),
                'longitude': ee.List(ee.List(triple).get(0)).get(1),
                'latitude': ee.List(triple).get(1)
            }
        )
    )
  ")

  unpack_func <- reticulate::py$unpack_compressed_geom
  unpacked_list <- compressed_fc$toList(num_feats)$map(unpack_func)$flatten()

  return(ee$FeatureCollection(unpacked_list))
}

#' GEE Classifier Methods Registry
#' @keywords internal
GEE_CLASSIFIER_METHODS <- list(
  rf = "smileRandomForest",
  gbt = "smileGradientTreeBoost",
  maxent = "amnhMaxent",
  svm = "libsvm"
)

#' GEE Reducer Methods Registry
#' @keywords internal
GEE_REDUCER_METHODS <- list(
  centroid = list(reducer = "mean", encoding = "presence"),
  ridge = list(reducer = "ridgeRegression", encoding = "presence")
)

#' Train GEE Model
#'
#' @param sampled_fc Sampled FeatureCollection (from get_embeddings_at_fc)
#' @param method Method name
#' @param params List of hyperparameters
#' @param class_property Column name for labels
#' @keywords internal
train_gee_model <- function(sampled_fc, method, params = list(), class_property = "present") {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)

  is_classifier <- method %in% names(GEE_CLASSIFIER_METHODS)
  is_reducer <- method %in% names(GEE_REDUCER_METHODS)

  if (!is_classifier && !is_reducer) stop("Unsupported method: ", method)

  LABEL_COL <- "label"

  # 1. Label Encoding
  if (is_reducer) {
    encoding <- GEE_REDUCER_METHODS[[method]]$encoding
    if (encoding == "binary") {
      # +1 / -1
      sampled_fc <- sampled_fc$map(function(f) {
        val <- ee$Number(f$get(class_property))
        lbl <- ee$Algorithms$If(val$eq(1), 1.0, -1.0)
        f$set(LABEL_COL, lbl)
      })
    } else {
      # 1 / 0
      sampled_fc <- sampled_fc$map(function(f) {
        f$set(LABEL_COL, ee$Number(f$get(class_property)))
      })
    }
  } else {
    # Classifier: 0/1 integer
    sampled_fc <- sampled_fc$map(function(f) {
      f$set(LABEL_COL, ee$Number(f$get(class_property))$toInt())
    })
  }

  # 2. Training
  if (is_classifier) {
    # Class balancing (1:1)
    pos_fc <- sampled_fc$filter(ee$Filter$eq(LABEL_COL, 1L))
    neg_fc <- sampled_fc$filter(ee$Filter$eq(LABEL_COL, 0L))
    pos_count <- pos_fc$size()
    balanced_fc <- pos_fc$merge(neg_fc$randomColumn()$sort("random")$limit(pos_count))

    gee_method <- GEE_CLASSIFIER_METHODS[[method]]
    # Build classifier with params
    clf_factory <- ee$Classifier[[gee_method]]

    # Map R params to GEE params based on method
    filtered_params <- list()
    if (gee_method == "smileRandomForest") {
      filtered_params <- list(
        numberOfTrees = as.integer(params$numberOfTrees),
        variablesPerSplit = if (!is.null(params$variablesPerSplit)) as.integer(params$variablesPerSplit) else NULL,
        minLeafPopulation = if (!is.null(params$minLeafPopulation)) as.integer(params$minLeafPopulation) else NULL,
        bagFraction = params$bagFraction,
        maxNodes = if (!is.null(params$maxNodes)) as.integer(params$maxNodes) else NULL
      )
    } else if (gee_method == "smileGradientTreeBoost") {
      filtered_params <- list(
        numberOfTrees = as.integer(params$numberOfTrees),
        shrinkage = params$shrinkage,
        maxNodes = if (!is.null(params$maxNodes)) as.integer(params$maxNodes) else NULL
      )
    } else {
      # For MaxEnt, SVM, etc., only pass what was in the user-provided 'options'
      # which we stored inside 'params' but we need to isolate them.
      # To keep it simple, if it's not one of the core tree params or lambda, pass it.
      core_params <- c("numberOfTrees", "minLeafPopulation", "bagFraction", "shrinkage", "maxNodes", "variablesPerSplit", "lambda_", "polynomial")
      filtered_params <- params[setdiff(names(params), core_params)]
    }

    # Remove NULLs
    filtered_params <- filtered_params[!sapply(filtered_params, is.null)]

    clf <- do.call(clf_factory, filtered_params)
    clf <- clf$setOutputMode("PROBABILITY")

    trained_model <- clf$train(
      features = balanced_fc,
      classProperty = LABEL_COL,
      inputProperties = emb_cols
    )

    # Retrieve model explanation (n_trees, variable importance, etc) for transparency
    explanation <- tryCatch(trained_model$explain()$getInfo(), error = function(e) list())

    return(list(
      trained = trained_model,
      is_classifier = TRUE,
      method = method,
      n_samples_total = as.integer(sampled_fc$size()$getInfo()),
      n_samples_balanced = as.integer(balanced_fc$size()$getInfo()),
      gee_explain = explanation
    ))
  } else {
    # Reducer logic (Ridge, Linear, Centroid)
    reducer_info <- GEE_REDUCER_METHODS[[method]]
    reducer_name <- reducer_info$reducer

    if (reducer_name %in% c("ridgeRegression", "linearRegression", "robustLinearRegression")) {
      lambda_val <- if (!is.null(params$lambda_)) params$lambda_ else 0.0
      poly_degree <- if (!is.null(params$polynomial)) as.integer(params$polynomial) else 1L

      training_fc <- sampled_fc
      current_emb_cols <- emb_cols

      if (poly_degree > 1) {
        # Quadratic: Linear + Squared terms
        training_fc <- training_fc$map(function(f) {
          arr <- f$toArray(emb_cols)
          sq <- arr$multiply(arr)
          dict <- ee$Dictionary$fromLists(paste0(emb_cols, "_2"), sq$toList())
          return(f$set(dict))
        })
        current_emb_cols <- c(emb_cols, paste0(emb_cols, "_2"))
      }

      num_x <- as.integer(length(current_emb_cols))
      reducer <- ee$Reducer$ridgeRegression(numX = num_x, numY = 1L, lambda = lambda_val)
      selectors <- c(current_emb_cols, LABEL_COL)

      result <- training_fc$reduceColumns(reducer = reducer, selectors = selectors)$getInfo()
      coefs <- as.numeric(result$coefficients)

      # GEE convention for these reducers is typically Intercept first
      intercept <- coefs[1]
      weights <- coefs[2:length(coefs)]

      return(list(
        weights = weights,
        intercept = intercept,
        is_classifier = FALSE,
        method = method,
        polynomial = poly_degree,
        n_samples = as.integer(sampled_fc$size()$getInfo())
      ))
    } else if (reducer_name == "mean") {
      # Centroid
      presence_fc <- sampled_fc$filter(ee$Filter$eq(LABEL_COL, 1.0))
      res <- presence_fc$reduceColumns(
        reducer = ee$Reducer$mean()$`repeat`(64L),
        selectors = emb_cols
      )$getInfo()
      weights <- as.numeric(res$mean)
      return(list(
        weights = weights,
        intercept = 0.0,
        is_classifier = FALSE,
        method = method,
        n_samples = as.integer(sampled_fc$size()$getInfo()),
        n_samples_presence = as.integer(presence_fc$size()$getInfo())
      ))
    }
  }
}

#' Assign Spatial Folds
#'
#' @param df Data frame with longitude, latitude
#' @param n_folds Number of folds
#' @keywords internal
assign_spatial_folds <- function(df, n_folds = 10) {
  ee <- reticulate::import("ee")
  if (nrow(df) <= n_folds) {
    df$fold <- (seq_len(nrow(df)) - 1) %% n_folds
    return(df)
  }

  # 1. Upload coordinates to GEE
  df$tmp_id_spatial <- seq_len(nrow(df)) - 1

  features <- lapply(seq_len(nrow(df)), function(i) {
    geom <- ee$Geometry$Point(c(as.numeric(df$longitude[i]), as.numeric(df$latitude[i])))
    ee$Feature(geom, list(
      tmp_id = as.integer(df$tmp_id_spatial[i]),
      lat = as.numeric(df$latitude[i]),
      lon = as.numeric(df$longitude[i])
    ))
  })

  fc <- ee$FeatureCollection(features)

  # 2. KMeans on GEE
  clusterer <- ee$Clusterer$wekaKMeans(as.integer(n_folds))$train(fc, list("lat", "lon"))
  clustered <- fc$cluster(clusterer, "fold")

  # 3. Retrieve
  cl_res <- clustered$reduceColumns(ee$Reducer$toList(2L), list("tmp_id", "fold"))$getInfo()

  # 4. Join
  fold_list <- cl_res$list
  # fold_list is a list of lists [[id, fold], ...]
  ids <- sapply(fold_list, `[[`, 1)
  folds <- sapply(fold_list, `[[`, 2)
  fold_map <- setNames(folds, as.character(ids))

  df$fold <- as.integer(fold_map[as.character(df$tmp_id_spatial)])
  df$tmp_id_spatial <- NULL

  return(df)
}

#' Generate Background Embeddings
#'
#' @param data Training data frame containing longitude, latitude, and year
#' @param scale Resolution in meters
#' @param count Total number of points to generate
#' @keywords internal
get_background_embeddings_gee <- function(data, scale, count) {
  ee <- reticulate::import("ee")

  bbox <- c(
    min(data$longitude, na.rm = TRUE), min(data$latitude, na.rm = TRUE),
    max(data$longitude, na.rm = TRUE), max(data$latitude, na.rm = TRUE)
  )

  # if there's only one point, buffer it to create a valid bbox
  if (bbox[1] == bbox[3] && bbox[2] == bbox[4]) {
    aoi_geom <- ee$Geometry$Point(c(bbox[1], bbox[2]))$buffer(1000)
    timestamp_message(sprintf("  -> Using buffered 1km point (single point dataset): [%.4f, %.4f]", bbox[1], bbox[2]))
  } else {
    aoi_geom <- ee$Geometry$Rectangle(bbox)
    timestamp_message(sprintf("  -> Using training data bounding box: [%.4f, %.4f, %.4f, %.4f]", bbox[1], bbox[2], bbox[3], bbox[4]))
  }

  year_counts <- table(data$year)
  total_train <- sum(year_counts)

  # Safe batch limit per single GEE sample() call based on resolution
  batch_limit <- if (scale <= 10) 1000L else if (scale <= 160) 5000L else 10000L

  timestamp_message("  -> Generating background points with embeddings (100% server-side)...")

  bg_fcs <- list()
  total_requested <- 0
  total_generated <- 0

  for (yr_str in names(year_counts)) {
    yr <- as.integer(yr_str)
    n_yr <- round(count * (year_counts[[yr_str]] / total_train))
    if (n_yr == 0) next

    img <- get_embedding_image(yr, scale)

    # Process one year at a time
    remaining <- n_yr
    batch_seed <- 0L
    yr_generated <- 0
    yr_fcs <- list()
    
    # Try up to 5x the requested count to find enough land points
    max_total_requests <- n_yr * 5
    attempts_made <- 0
    
    while (remaining > 0 && attempts_made < max_total_requests) {
      n_batch <- min(remaining * 2, batch_limit) # Request 2x remaining to fill gaps faster
      attempts_made <- attempts_made + n_batch
      
      sampled <- img$sample(
        region = aoi_geom,
        scale = scale,
        numPixels = as.integer(n_batch),
        geometries = FALSE,
        tileScale = 16L,
        seed = as.integer(yr * 100L + batch_seed)
      )$filter(ee$Filter$notNull(list("A00")))$
        map(function(f) f$set("present", 0L)$set("year", as.integer(yr)))

      batch_count <- as.integer(sampled$size()$getInfo())
      if (batch_count > 0) {
        # Take only what we need to avoid exceeding the target
        take_now <- min(batch_count, remaining)
        if (take_now < batch_count) {
            sampled <- sampled$limit(as.integer(take_now))
        }
        
        yr_fcs <- c(yr_fcs, list(sampled))
        yr_generated <- yr_generated + take_now
        remaining <- remaining - take_now
      }
      
      batch_seed <- batch_seed + 1L
      if (batch_count == 0 && attempts_made < max_total_requests) {
        # If we got 0, the AOI might be mostly water. 
        # We'll try one more large batch before giving up on this year.
        n_batch <- batch_limit 
      }
    }
    
    if (length(yr_fcs) > 0) {
        total_requested <- total_requested + n_yr
        total_generated <- total_generated + yr_generated
        timestamp_message(sprintf("     Year %d: Target %d, Final %d (Attempts: %d points sampled)", yr, n_yr, yr_generated, attempts_made))
        
        # Merge all chunks for this year
        yr_merged <- yr_fcs[[1]]
        if (length(yr_fcs) > 1) {
            for (j in 2:length(yr_fcs)) yr_merged <- yr_merged$merge(yr_fcs[[j]])
        }
        bg_fcs <- c(bg_fcs, list(yr_merged))
    }
  }

  # Final Merge
  bg_fc <- bg_fcs[[1]]
  if (length(bg_fcs) > 1) {
    for (j in 2:length(bg_fcs)) {
      bg_fc <- bg_fc$merge(bg_fcs[[j]])
    }
  }

  timestamp_message(sprintf("  -> Total Background Points Finalized: %d\n", total_generated))

  return(bg_fc)
}

#' Predict GEE Map
#'
#' @param model_res Result from train_gee_model
#' @param img Alpha Earth mosaic
#' @keywords internal
predict_gee_map <- function(model_res, img) {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)

  if (model_res$is_classifier) {
    score_col <- if (model_res$method == "maxent") "probability" else "classification"
    prediction <- img$classify(model_res$trained)$select(score_col)

    if (model_res$method == "svm") {
      prediction <- ee$Image(1.0)$subtract(prediction)
    }

    return(prediction$rename("similarity"))
  } else {
    # Reducer dot product
    weights <- model_res$weights
    intercept <- model_res$intercept
    poly_degree <- if (!is.null(model_res$polynomial)) model_res$polynomial else 1L

    input_img <- img
    current_cols <- emb_cols
    if (poly_degree > 1) {
      img_sq <- img$multiply(img)$rename(paste0(emb_cols, "_2"))
      input_img <- img$addBands(img_sq)
      current_cols <- c(emb_cols, paste0(emb_cols, "_2"))
    }

    weights_ee <- ee$Image$constant(as.list(weights))$rename(current_cols)
    prediction <- input_img$multiply(weights_ee)$reduce(ee$Reducer$sum())$add(intercept)

    return(prediction$rename("similarity"))
  }
}

#' Predict Scores for Multiple Models on GEE
#'
#' Consolidates scoring to minimize GEE map passes.
#'
#' @param fc GEE FeatureCollection with embeddings
#' @param models_list List of model results from AlphaSDM
#' @keywords internal
predict_all_models_gee <- function(fc, models_list) {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)

  # 1. Identify Classifiers vs Reducers
  methods <- names(models_list)
  classifiers <- methods[sapply(models_list, function(m) m$is_classifier)]
  reducers <- methods[!sapply(models_list, function(m) m$is_classifier)]

  scored_fc <- fc

  # 2. Chain Classifiers (Native and fast)
  for (m in classifiers) {
    model_res <- models_list[[m]]
    score_col <- if (model_res$method == "maxent") "probability" else "classification"

    # Run native classification
    scored_fc <- scored_fc$classify(model_res$trained)

    # Handle SVM inversion or MaxEnt naming
    target_col <- paste0("pred_", m)
    if (model_res$method == "svm") {
      scored_fc <- scored_fc$map(function(f) f$set(target_col, ee$Number(1.0)$subtract(f$get(score_col))))
    } else {
      scored_fc <- scored_fc$map(function(f) f$set(target_col, f$get(score_col)))
    }
  }

  # 3. Consolidate Reducers into ONE map pass (Math expression map)
  if (length(reducers) > 0) {
    scored_fc <- scored_fc$map(function(f) {
      emb <- ee$Array(f$toArray(emb_cols))

      # Prepare squared terms if any model needs them
      has_poly <- any(sapply(models_list[reducers], function(m) isTRUE(m$polynomial > 1)))
      emb_aug <- emb
      if (has_poly) {
        sq <- emb$multiply(emb)
        emb_aug <- ee$Array$cat(list(emb, sq), 0L)
      }

      # Scored features container
      f_out <- f
      for (m in reducers) {
        model_res <- models_list[[m]]
        poly_degree <- if (!is.null(model_res$polynomial)) model_res$polynomial else 1L

        # Select active part of augmented array
        current_emb <- if (poly_degree > 1) emb_aug else emb
        weights_ee <- ee$Array(as.list(model_res$weights))
        intercept <- model_res$intercept

        # Dot product
        score <- current_emb$multiply(weights_ee)$reduce(ee$Reducer$sum(), list(0L))$get(list(0L))$add(intercept)
        f_out <- f_out$set(paste0("pred_", m), score)
      }
      return(f_out)
    })
  }

  return(scored_fc)
}

#' Predict Scores for Single Model on GEE
#'
#' @param fc GEE FeatureCollection with embeddings
#' @param model_res Model result
#' @keywords internal
predict_points_gee <- function(fc, model_res) {
  # Wrapper for backward compatibility
  models_list <- list()
  models_list[[model_res$method]] <- model_res
  res_fc <- predict_all_models_gee(fc, models_list)

  # Map back to standard "score" property for caller
  predicted_col <- paste0("pred_", model_res$method)
  return(res_fc$map(function(f) f$set("score", f$get(predicted_col))))
}

#' Internal Training Pipeline
#'
#' @keywords internal
train_AlphaSDM_internal <- function(data, methods, aoi_geom = NULL, scale = 10, aoi_year = NULL,
                                    count = NULL, training_params = list()) {
  # 1. Background Logic
  needs_bg <- any(methods %in% c("centroid", "ridge", "rf", "gbt", "maxent", "svm")) &&
    !("present" %in% names(data) && any(data$present == 0))

  sampled_fc <- NULL
  if (needs_bg) {
    timestamp_message("\n--- Generating Background Points (GEE-side) ---")
    n_bg <- if (!is.null(count)) count else nrow(data)
    bg_fc <- get_background_embeddings_gee(data, scale, n_bg)

    # Prepare presence points
    prep_pres <- prep_training_data_gee(data, class_property = "present", scale = scale)
    pres_fc <- prep_pres$fc

    sampled_fc <- pres_fc$merge(bg_fc)
  } else {
    prep <- prep_training_data_gee(data, class_property = "present", scale = scale)
    sampled_fc <- prep$fc
  }

  # 2. Training Loop
  timestamp_message(sprintf("--- Training %d Method(s) via R-GEE ---", length(methods)))

  models <- list()
  meta_results <- list()

  for (m in methods) {
    timestamp_message(sprintf("  Training %s...", m))

    start_t <- Sys.time()
    models[[m]] <- train_gee_model(sampled_fc, m, params = training_params)
    end_t <- Sys.time()

    # Store metadata
    meta <- models[[m]]
    meta$trained <- NULL # Remove heavy GEE object
    meta$training_seconds <- as.numeric(difftime(end_t, start_t, units = "secs"))
    meta_results[[m]] <- meta
  }

  return(list(models = models, metadata = meta_results))
}
