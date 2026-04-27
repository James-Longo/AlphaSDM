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
    years <- retry_curl_download(fc$aggregate_array("year")$distinct()$getInfo())
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

#' Upload Points to GEE efficiently
#' @param df Data frame with longitude, latitude, and optional present column.
#' @return ee$FeatureCollection
upload_points_to_gee <- function(df) {
  ee <- reticulate::import("ee")
  
  # Prepare features with properties dynamically
  features <- lapply(seq_len(nrow(df)), function(i) {
    props <- as.list(df[i, , drop = FALSE])
    # ee$Feature expects a list of primitives (scalars)
    props_clean <- lapply(props, function(x) {
      val <- if (is.numeric(x)) as.numeric(x) else if (is.integer(x)) as.integer(x) else x
      if (length(val) > 0) val[1] else NULL
    })
    
    ee$Feature(
      ee$Geometry$Point(c(as.numeric(props$longitude), as.numeric(props$latitude))),
      props_clean
    )
  })
  
  return(ee$FeatureCollection(features))
}

#' Prepare Training Data for GEE
prep_training_data_gee <- function(df, class_property = "present", scale = 10) {
  ee <- reticulate::import("ee")

  # 1. Clean data
  cols_to_keep <- c("longitude", "latitude", "year")
  if (!is.null(class_property)) cols_to_keep <- c(cols_to_keep, class_property)

  df_clean <- df[complete.cases(df[, cols_to_keep]), ]
  if (nrow(df_clean) == 0) stop("No valid training data remaining after dropping NAs.")

  # 2. Upload to GEE
  timestamp_message("\n--- Uploading User Data ---")
  timestamp_message(sprintf("  -> Transferring %d provided coordinates to GEE...", nrow(df_clean)))

  # Use the specialized uploader for training data
  upload_fc <- upload_points_to_gee(df_clean)

  # 3. Sample embeddings
  sample_props <- c("year", class_property)
  sampled_fc <- get_embeddings_at_fc(upload_fc, scale, properties = sample_props)

  # Log valid points remaining
  final_count <- as.numeric(retry_curl_download(sampled_fc$size()$getInfo()))
  dropped <- nrow(df_clean) - final_count
  if (dropped > 0) {
    timestamp_message(sprintf("  -> Discarded %d presence points due to missing embeddings (e.g., over water).", dropped))
  }
  timestamp_message(sprintf("  -> Final valid presence points: %d", final_count))

  return(list(
    fc = sampled_fc,
    df_clean = df_clean,
    class_property = class_property
  ))
}

#' GEE Classifier Methods Registry
GEE_CLASSIFIER_METHODS <- list(
  rf = "smileRandomForest",
  gbt = "smileGradientTreeBoost",
  maxent = "amnhMaxent",
  svm = "libsvm"
)

#' GEE Reducer Methods Registry
GEE_REDUCER_METHODS <- list(
  similarity = list(reducer = "mean", encoding = "presence")
)

#' Train GEE Model
train_gee_model <- function(sampled_fc, method, params = list(), class_property = "present") {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)

  is_classifier <- method %in% names(GEE_CLASSIFIER_METHODS)
  is_reducer <- method %in% names(GEE_REDUCER_METHODS)

  if (!is_classifier && !is_reducer) stop("Unsupported method: ", method)

  LABEL_COL <- "label"

  # 1. Label Encoding
  sampled_fc <- sampled_fc$map(function(f) {
    f$set(LABEL_COL, ee$Number(f$get(class_property))$toInt())
  })

  # 2. Training
  if (is_classifier) {
    # Class balancing (1:1)
    pos_fc <- sampled_fc$filter(ee$Filter$eq(LABEL_COL, 1L))
    neg_fc <- sampled_fc$filter(ee$Filter$eq(LABEL_COL, 0L))
    pos_count <- pos_fc$size()
    balanced_fc <- pos_fc$merge(neg_fc$randomColumn("random", 42L)$sort("random")$limit(pos_count))

    gee_method <- GEE_CLASSIFIER_METHODS[[method]]
    clf_factory <- ee$Classifier[[gee_method]]

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
      core_params <- c("numberOfTrees", "minLeafPopulation", "bagFraction", "shrinkage", "maxNodes", "variablesPerSplit", "lambda_", "polynomial", "batch_size")
      filtered_params <- params[setdiff(names(params), core_params)]
    }

    filtered_params <- filtered_params[!sapply(filtered_params, is.null)]
    clf <- do.call(clf_factory, filtered_params)
    clf <- clf$setOutputMode("PROBABILITY")

    trained_model <- clf$train(
      features = balanced_fc,
      classProperty = LABEL_COL,
      inputProperties = emb_cols
    )

    explanation <- tryCatch(retry_curl_download(trained_model$explain()$getInfo()), error = function(e) list())

    return(list(
      trained = trained_model,
      is_classifier = TRUE,
      method = method,
      n_samples_total = as.integer(sampled_fc$size()$getInfo()),
      n_samples_balanced = as.integer(balanced_fc$size()$getInfo()),
      gee_explain = explanation
    ))
  } else {
    # Reducer logic (Simplified Centroid)
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
      n_samples = as.numeric(sampled_fc$size()$getInfo())
    ))
  }
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
    # Simple Dot Product Similarity
    weights_ee <- ee$Image$constant(as.list(model_res$weights))$rename(emb_cols)
    prediction <- img$multiply(weights_ee)$reduce(ee$Reducer$sum())
    
    return(prediction$rename("similarity"))
  }
}

#' Predict Scores for Multiple Models on GEE
#'
#' @param fc GEE FeatureCollection with embeddings
#' @param models_list List of model results from AlphaSDM
#' @keywords internal
predict_all_models_gee <- function(fc, models_list) {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)

  methods <- names(models_list)
  classifiers <- methods[sapply(models_list, function(m) m$is_classifier)]
  reducers <- methods[!sapply(models_list, function(m) m$is_classifier)]

  scored_fc <- fc

  # 1. Chain Classifiers
  for (m in classifiers) {
    model_res <- models_list[[m]]
    score_col <- if (model_res$method == "maxent") "probability" else "classification"
    scored_fc <- scored_fc$classify(model_res$trained)
    
    target_col <- paste0("pred_", m)
    if (model_res$method == "svm") {
      scored_fc <- scored_fc$map(function(f) f$set(target_col, ee$Number(1.0)$subtract(f$get(score_col))))
    } else {
      scored_fc <- scored_fc$map(function(f) f$set(target_col, f$get(score_col)))
    }
  }

  # 2. Handle Reducers (Similarity)
  if (length(reducers) > 0) {
    for (m in reducers) {
        model_res <- models_list[[m]]
        centroid_ee <- ee$Array(as.list(model_res$weights))
        
        target_col <- paste0("pred_", m)
        scored_fc <- scored_fc$map(function(f) {
            point_ee <- ee$Array(f$toArray(emb_cols))
            # Linear dot product
            score <- point_ee$multiply(centroid_ee)$reduce(ee$Reducer$sum(), list(0L))$get(list(0L))
            return(f$set(target_col, score))
        })
    }
  }

  return(scored_fc)
}

#' Assign Spatial Folds
assign_spatial_folds <- function(df, n_folds = 10) {
  ee <- reticulate::import("ee")
  if (nrow(df) <= n_folds) {
    df$fold <- (seq_len(nrow(df)) - 1) %% n_folds
    return(df)
  }

  df$tmp_id <- seq_len(nrow(df)) - 1
  if (!"year" %in% names(df)) df$year <- 2020
  fc <- upload_points_to_gee(df)

  clusterer <- ee$Clusterer$wekaKMeans(as.integer(n_folds))$train(fc, list("latitude", "longitude"))
  clustered <- fc$cluster(clusterer, "fold")

  cl_res <- retry_curl_download(clustered$reduceColumns(ee$Reducer$toList(2L), list("tmp_id", "fold"))$getInfo())

  fold_list <- cl_res$list
  ids <- sapply(fold_list, `[[`, 1)
  folds <- sapply(fold_list, `[[`, 2)
  fold_map <- setNames(folds, as.character(ids))

  df$fold <- as.integer(fold_map[as.character(df$tmp_id)])
  df$tmp_id <- NULL

  return(df)
}

#' Get Feature Names
#' @keywords internal
get_feature_names <- function(fc) {
  return(sprintf("A%02d", 0:63))
}
