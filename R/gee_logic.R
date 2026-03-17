#' Get Alpha Earth Embedding Image
#' 
#' @param year Year (2017-2025)
#' @param scale Resolution in meters
#' @keywords internal
get_embedding_image <- function(year, scale) {
  ee <- reticulate::import("ee")
  asset_path <- "GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL"
  emb_cols <- sprintf("A%02d", 0:63)
  
  img <- ee$ImageCollection(asset_path)$
    filter(ee$Filter$calendarRange(as.integer(year), as.integer(year), "year"))$
    mosaic()$
    select(emb_cols)
    
  if (scale > 10) {
    # Calculate exact maxPixels needed with a 20% safety buffer
    # This is scientifically accurate mean reduction (average of pixels in footprint)
    pixel_ratio <- (scale / 10)^2
    safe_max_pixels <- ceiling(pixel_ratio * 1.5)
    
    img <- img$setDefaultProjection("EPSG:3857", NULL, 10)$reduceResolution(
      reducer = ee$Reducer$mean(),
      maxPixels = safe_max_pixels
    )$reproject(
      crs = "EPSG:3857",
      scale = scale
    )
  }
  
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
get_embeddings_at_fc <- function(fc, scale, properties = NULL, geometries = FALSE, years = NULL) {
  ee <- reticulate::import("ee")
  
  if (is.null(years)) {
    #aggregate_array().getInfo() can be expensive for large FCs
    years <- fc$aggregate_array("year")$distinct()$getInfo()
  }
  t0 <- log_time(sprintf("Sampling embeddings at %d years", length(years)))
  sampled_fcs <- list()
  for (yr in years) {
    yr_start <- log_time(sprintf("  Processing year %d", yr))
    yr_fc <- fc$filter(ee$Filter$eq("year", as.integer(yr)))
    yr_img <- get_embedding_image(yr, scale)
    
    sampled <- yr_img$sampleRegions(
      collection = yr_fc,
      properties = as.list(properties),
      scale = scale,
      geometries = geometries,
      tileScale = 16L
    )$filter(ee$Filter$notNull(as.list("A00")))
    
    sampled_fcs <- c(sampled_fcs, list(sampled))
    log_time(sprintf("  Year %d", yr), yr_start)
  }
  
  res <- ee$FeatureCollection(sampled_fcs)$flatten()
  log_time("Sampling complete", t0)
  return(res)
}

#' Prepare Training Data for GEE
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
  
  # 2. Upload to GEE as MultiPoint features for efficiency
  # We group by year and class to minimize feature count
  group_cols <- "year"
  if (!is.null(class_property)) group_cols <- c(group_cols, class_property)
  
  # For R, we can use interaction for multi-column grouping
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
    
    coords_mat <- as.matrix(grp[, c("longitude", "latitude")])
    
    props <- list(year = as.integer(yr))
    if (!is.null(class_property)) props[[class_property]] <- cls_val
    
    # Chunk MultiPoints if they get too large
    chunk_size <- 5000
    for (i in seq(1, nrow(coords_mat), by = chunk_size)) {
      end <- min(i + chunk_size - 1, nrow(coords_mat))
      sub_coords <- coords_mat[i:end, , drop = FALSE]
      
      # Convert matrix to list of lists for GEE Python compatibility
      coord_list <- lapply(seq_len(nrow(sub_coords)), function(j) as.list(unname(sub_coords[j, ])))
      
      geom <- ee$Geometry$MultiPoint(coord_list)
      gee_features <- c(gee_features, list(ee$Feature(geom, props)))
    }
  }
  
  t_agg <- log_time("Aggregating points for GEE")
  gee_features_list <- ee$FeatureCollection(gee_features)
  log_time("Aggregation complete", t_agg)
  
  # 3. Sample embeddings
  t_samp_inner <- log_time(sprintf("Sampling Alpha Earth (%dm footprint mean)", scale))
  sampled_fc <- get_embeddings_at_fc(gee_features_list, scale, properties = group_cols)
  log_time("Sampling complete", t_samp_inner)
  
  return(list(
    fc = sampled_fc,
    df_clean = df_clean,
    class_property = class_property
  ))
}

#' GEE Classifier Methods Registry
#' @keywords internal
GEE_CLASSIFIER_METHODS <- list(
  rf = "smileRandomForest",
  gbt = "smileGradientTreeBoost",
  cart = "smileCart",
  maxent = "amnhMaxent",
  svm = "libsvm"
)

#' GEE Reducer Methods Registry
#' @keywords internal
GEE_REDUCER_METHODS <- list(
  centroid = list(reducer = "mean", encoding = "presence"),
  ridge = list(reducer = "ridgeRegression", encoding = "binary"),
  linear = list(reducer = "linearRegression", encoding = "binary"),
  robust_linear = list(reducer = "robustLinearRegression", encoding = "binary")
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
    } else if (gee_method == "smileCart") {
      filtered_params <- list(
        maxNodes = if (!is.null(params$maxNodes)) as.integer(params$maxNodes) else NULL,
        minLeafPopulation = if (!is.null(params$minLeafPopulation)) as.integer(params$minLeafPopulation) else NULL
      )
    } else {
      # For others, try to pass whatever is explicitly requested in 'options' 
      # or just stay minimal.
      filtered_params <- list()
    }
    
    # Remove NULLs
    filtered_params <- filtered_params[!sapply(filtered_params, is.null)]
    
    clf <- do.call(clf_factory, filtered_params)
    clf <- clf$setOutputMode("PROBABILITY")
    
    t_train <- log_time(sprintf("Training GEE model: %s", method))
    trained_model <- clf$train(
      features = balanced_fc,
      classProperty = LABEL_COL,
      inputProperties = emb_cols
    )
    log_time("Training complete", t_train)
    
    return(list(trained = trained_model, is_classifier = TRUE, method = method))
    
  } else {
    # Reducer logic (Ridge, Linear, Centroid)
    reducer_info <- GEE_REDUCER_METHODS[[method]]
    reducer_name <- reducer_info$reducer
    
    if (reducer_name %in% c("ridgeRegression", "linearRegression", "robustLinearRegression")) {
      # Case weight balancing handled in R before upload? 
      # Actually Python handled it per-fold or per-run.
      # For now, let's assume balanced or weighted.
      # (Skipping centering logic for now to keep it simple, can add later if accuracy differs)
      
      if (reducer_name == "ridgeRegression") {
        lambda <- if (!is.null(params$lambda_)) params$lambda_ else 0.1
        reducer <- ee$Reducer$ridgeRegression(numX = 64L, numY = 1L, lambda = lambda)
      } else if (reducer_name == "robustLinearRegression") {
        beta <- if (!is.null(params$beta)) params$beta else NULL
        reducer <- ee$Reducer$robustLinearRegression(numX = 65L, numY = 1L, beta = beta)
        # Note: robust needs constant column
      } else {
        reducer <- ee$Reducer$linearRegression(numX = 65L, numY = 1L)
      }
      
      # For linear/robust, need constant column
      if (reducer_name != "ridgeRegression") {
        sampled_fc <- sampled_fc$map(function(f) f$set("constant", 1.0))
        selectors <- c("constant", emb_cols, LABEL_COL)
      } else {
        selectors <- c(emb_cols, LABEL_COL)
      }
      
      result <- sampled_fc$reduceColumns(reducer = reducer, selectors = selectors)$getInfo()
      coefs <- as.numeric(result$coefficients)
      
      if (reducer_name == "ridgeRegression") {
        # Ridge GEE returns [intercept, w0, w1, ...]
        intercept <- coefs[1]
        weights <- coefs[2:65]
      } else {
        intercept <- coefs[1]
        weights <- coefs[2:65]
      }
      
      return(list(weights = weights, intercept = intercept, is_classifier = FALSE, method = method))
      
    } else if (reducer_name == "mean") {
      # Centroid
      presence_fc <- sampled_fc$filter(ee$Filter$eq(LABEL_COL, 1.0))
      
      # We calculate the mean using a chunked-sum approach on GEE
      # to avoid memory limits and ENSURE no raw embeddings are pulled to R.
      t_cent <- log_time("Calculating Centroid weights (Server-side Aggregate)")
      
      n_presence <- as.numeric(presence_fc$size()$getInfo())
      if (n_presence == 0) stop("Centroid training failed: No presence points found.")
      
      chunk_size_agg <- 5000
      full_list <- presence_fc$toList(n_presence)
      total_sum <- numeric(64)
      
      for (i in seq(0, n_presence - 1, by = chunk_size_agg)) {
          end_idx <- min(i + chunk_size_agg, n_presence)
          chunk_fc <- ee$FeatureCollection(full_list$slice(i, end_idx))
          
          # Compute sum of embeddings for this chunk on GEE
          chunk_sum <- chunk_fc$reduceColumns(
              reducer = ee$Reducer$sum()$`repeat`(64L),
              selectors = as.list(emb_cols)
          )$getInfo()
          
          total_sum <- total_sum + as.numeric(chunk_sum$sum)
      }
      
      weights <- total_sum / n_presence
      log_time("Centroid calculation", t_cent)
      
      return(list(weights = as.numeric(weights), intercept = 0.0, is_classifier = FALSE, method = method, trained = NULL))
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
    ee$Feature(geom, list(tmp_id = as.integer(df$tmp_id_spatial[i]), 
                           lat = as.numeric(df$latitude[i]), 
                           lon = as.numeric(df$longitude[i])))
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
#' @param aoi_geom GEE Geometry
#' @param year Year for Alpha Earth mosaic
#' @param scale Resolution in meters
#' @param count Number of points
#' @keywords internal
get_background_embeddings_gee <- function(aoi_geom, year, scale, count) {
  img <- get_embedding_image(year, scale)
  
  # Sample background points from the image within AOI
  sampled <- img$sample(
    region = aoi_geom,
    scale = scale,
    numPixels = as.integer(count),
    geometries = TRUE,
    tileScale = 16L
  )$map(function(f) f$set("present", 0L)$set("year", as.integer(year)))
  
  return(sampled)
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
    
    weights_ee <- ee$Image$constant(as.list(weights))$rename(emb_cols)
    prediction <- img$multiply(weights_ee)$reduce(ee$Reducer$sum())$add(intercept)
    
    return(prediction$rename("similarity"))
  }
}

#' Predict at coordinates entirely on GEE
#' @keywords internal
predict_at_coords_gee <- function(coords_df, models, methods, scale, gee_project = NULL) {
    ee <- reticulate::import("ee")
    t_tot <- log_time(sprintf("Predicting %d points entirely on GEE", nrow(coords_df)))
    
    # 1. Upload points
    t_prep <- log_time("Uploading points to GEE")
    all_features <- list()
    for (i in seq_len(nrow(coords_df))) {
        all_features[[i]] <- ee$Feature(
            ee$Geometry$Point(c(as.numeric(coords_df$longitude[i]), as.numeric(coords_df$latitude[i]))),
            list(year = as.integer(coords_df$year[i]), tmp_id = as.integer(i-1))
        )
    }
    upload_fc <- ee$FeatureCollection(all_features)
    log_time("Upload complete", t_prep)
    
    # 2. Sample embeddings (Still on GEE)
    t_samp <- log_time("Sampling embeddings (GEE-side)")
    sampled_fc <- get_embeddings_at_fc(upload_fc, scale, properties = c("year", "tmp_id"))
    log_time("Sampling complete", t_samp)
    
    # 3. Apply models
    current_fc <- sampled_fc
    pred_cols <- c()
    emb_cols <- sprintf("A%02d", 0:63)
    
    for (m in methods) {
        m_start <- log_time(sprintf("Applying %s on GEE", m))
        model_res <- models[[m]]
        pred_col <- paste0("pred_", m)
        pred_cols <- c(pred_cols, pred_col)
        
        if (model_res$is_classifier) {
            # Correct GEE usage: classify the entire collection
            current_fc <- current_fc$classify(model_res$trained, pred_col)
            
            # Specific fix for SVM which returns inverted in GEE
            if (m == "svm") {
                current_fc <- current_fc$map(function(f) {
                   f$set(pred_col, ee$Number(1.0)$subtract(f$getNumber(pred_col)))
                })
            }
        } else {
            # Dot product for reducers (Linear, Ridge, Centroid)
            weights <- as.numeric(model_res$weights)
            intercept <- as.numeric(model_res$intercept)
            
            # Map over FC to do dot product
            current_fc <- current_fc$map(function(f) {
                # Use array math for efficiency
                vals <- f$toArray(as.list(emb_cols))
                score <- vals$dotProduct(ee$Array(as.list(weights)))$add(intercept)
                f$set(pred_col, score)
            })
        }
        log_time(sprintf("Model %s complete", m), m_start)
    }
    
    # 4. Ensemble
    if (length(methods) > 1) {
        t_ens <- log_time("Generating Ensemble on GEE")
        current_fc <- current_fc$map(function(f) {
            scores <- f$toArray(as.list(pred_cols))
            f$set("pred_ensemble", scores$reduce(ee$Reducer$mean(), list(0L))$get(list(0L)))
        })
        pred_cols <- c(pred_cols, "pred_ensemble")
        log_time("Ensemble complete", t_ens)
    }
    
    # 5. Retrieve ONLY scores (Scientific integrity: preserve all points)
    t_ret <- log_time("Retrieving scores only (joining records)")
    
    # Strip heavy embeddings before retrieval to stay under GEE limits
    current_fc <- current_fc$select(as.list(c(pred_cols, "tmp_id")))
    
    res_list <- retrieve_gee_properties_chunked(current_fc, c(pred_cols, "tmp_id"))
    log_time("Retrieval complete", t_ret)
    
    if (length(res_list) == 0) {
        # Fallback to empty shell with NAs
        final_df <- data.frame(tmp_id = 0:(nrow(coords_df)-1))
        for (col in pred_cols) final_df[[col]] <- NA
    } else {
        res_df <- as.data.frame(do.call(rbind, res_list))
        colnames(res_df) <- c(pred_cols, "tmp_id")
        res_df$tmp_id <- as.numeric(res_df$tmp_id)
        
        # Left merge to ensure original size and order (NAs for masked points)
        final_df <- merge(data.frame(tmp_id = 0:(nrow(coords_df)-1)), res_df, by = "tmp_id", all.x = TRUE)
        final_df <- final_df[order(final_df$tmp_id), ]
    }
    
    final_df$tmp_id <- NULL
    final_df$longitude <- coords_df$longitude
    final_df$latitude <- coords_df$latitude
    
    log_time("Predicting complete", t_tot)
    return(final_df)
}
