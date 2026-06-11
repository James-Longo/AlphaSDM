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
      sdm_warn(sprintf("GEE rate limit or timeout — retrying in %.0fs (attempt %d of %d)",
                       delay, i, max_retries), indent = 1L)
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

#' Download Background Embeddings via Independent-Chunk sampleRegions + computeFeatures
#'
#' Generates background points in independent small batches (different seeds),
#' sampling embeddings for each via ee.data.computeFeatures() pagination.
#' Each batch is a fresh, isolated server-side expression — no large FC to count
#' or slice, no accumulated lazy graph. Repeats until n valid points collected.
#' Returns an R data frame with longitude, latitude, and A00-A63 columns.
#' @keywords internal
get_embeddings_as_df <- function(region, img, n, scale, seed = 42L,
                                 chunk_size = 5000L, page_size = 500L) {
  ee      <- reticulate::import("ee")
  ee_data <- reticulate::import("ee.data")
  emb_cols <- sprintf("A%02d", 0:63)

  all_dfs   <- list()
  total     <- 0L
  s         <- as.integer(seed)
  max_iters <- ceiling(n / chunk_size) * 5L  # generous safety limit

  for (iter in seq_len(max_iters)) {
    if (total >= n) break

    chunk_pts <- ee$FeatureCollection$randomPoints(
      region, as.integer(chunk_size), seed = s
    )
    s <- s + 1L

    expr <- img$sampleRegions(
      collection = chunk_pts,
      scale      = as.integer(scale),
      tileScale  = 16L,
      geometries = TRUE
    )

    tok <- NULL
    repeat {
      params <- list(expression = expr, pageSize = as.integer(page_size))
      if (!is.null(tok) && nchar(tok) > 0L) params$pageToken <- tok

      r <- retry_curl_download(ee_data$computeFeatures(params))

      feats <- r[["features"]]
      if (!is.null(feats) && length(feats) > 0L) {
        mat <- do.call(rbind, lapply(feats, function(f) {
          coords <- f$geometry$coordinates
          if (is.null(coords)) return(NULL)
          vals <- lapply(emb_cols, function(e) f$properties[[e]])
          if (any(vapply(vals, is.null, logical(1L)))) return(NULL)
          c(as.numeric(coords[[1L]]), as.numeric(coords[[2L]]), as.numeric(vals))
        }))
        if (!is.null(mat) && nrow(mat) > 0L) {
          df           <- as.data.frame(mat, stringsAsFactors = FALSE)
          colnames(df) <- c("longitude", "latitude", emb_cols)
          df           <- df[!is.na(df$A00), , drop = FALSE]
          total        <- total + nrow(df)
          all_dfs      <- c(all_dfs, list(df))
        }
      }

      tok <- r[["nextPageToken"]]
      if (is.null(tok) || nchar(tok) == 0L) break
      if (total >= n) break
    }
  }

  if (length(all_dfs) == 0L) return(data.frame())
  result <- do.call(rbind, all_dfs)
  result[seq_len(min(n, nrow(result))), , drop = FALSE]
}

get_embeddings_at_fc <- function(fc, scale, properties = NULL, geometries = FALSE, years = NULL) {
  ee <- reticulate::import("ee")
  raw <- get_embeddings_at_fc_raw(fc, scale, properties, geometries, years)
  return(raw$filter(ee$Filter$notNull(as.list("A00"))))
}

#' Upload Points to GEE efficiently, chunking large DFs to stay under 10 MB
#' @param df Data frame with longitude, latitude, and optional columns.
#' @param chunk_size Max rows per GeoJSON payload (default 5000 ≈ 4 MB).
#' @return ee$FeatureCollection
upload_points_to_gee <- function(df, chunk_size = 5000L) {
  ee       <- reticulate::import("ee")
  json_mod <- reticulate::import("json")

  upload_chunk <- function(chunk_df) {
    coord_cols <- c("longitude", "latitude")
    prop_cols  <- setdiff(names(chunk_df), coord_cols)
    features   <- vector("list", nrow(chunk_df))
    for (i in seq_len(nrow(chunk_df))) {
      props <- setNames(
        lapply(prop_cols, function(col) {
          v <- chunk_df[[col]][[i]]
          if (is.integer(v)) as.integer(v) else if (is.numeric(v)) as.numeric(v) else v
        }),
        prop_cols
      )
      features[[i]] <- list(
        type       = "Feature",
        geometry   = list(type = "Point",
                          coordinates = list(chunk_df$longitude[[i]], chunk_df$latitude[[i]])),
        properties = props
      )
    }
    geojson_py <- json_mod$loads(
      as.character(jsonlite::toJSON(list(type = "FeatureCollection", features = features),
                                    auto_unbox = TRUE, digits = 10))
    )
    ee$FeatureCollection(geojson_py)
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

#' Generate Background Points Server-Side
#'
#' Generates pseudo-absence points entirely on GEE. Samples only band A00
#' for land validation — never downloads background coordinates to R.
#' Returns a GEE FeatureCollection (geometry + year + present=0).
#' @keywords internal
generate_background_fc_gee <- function(bbox, aoi_year, count, aoi_geom = NULL) {
  ee     <- reticulate::import("ee")
  region <- if (!is.null(aoi_geom)) aoi_geom else ee$Geometry$Rectangle(bbox)

  sdm_info(sprintf("Target: %d background points (server-side, no pre-check)", count), indent = 1L)

  # No A00 pre-check needed — get_embeddings_at_fc already filters notNull("A00"),
  # so ocean/invalid points are dropped during the main 64-band embedding sampling.
  # This collapses two sampleRegions calls into one.
  raw_fc <- ee$FeatureCollection$randomPoints(region, as.integer(count * 3L))
  bg_fc  <- raw_fc$map(function(f) f$set("year", as.integer(aoi_year), "present", 0L))
  return(bg_fc$limit(as.integer(count)))
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
  sdm_section("Uploading training data to Google Earth Engine")
  pb_up <- sdm_progress_start("Uploading and sampling")
  sdm_info(sprintf("Transferring %d coordinates ...", nrow(df_clean)))

  # Use the specialized uploader for training data
  upload_fc <- upload_points_to_gee(df_clean)

  # 3. Sample embeddings
  sample_props <- c("year", class_property)
  known_years <- as.list(unique(as.integer(df_clean$year)))
  sampled_fc  <- get_embeddings_at_fc(upload_fc, scale, properties = sample_props, years = known_years)

  # Log valid points remaining (coverage filter)
  after_coverage <- as.numeric(retry_curl_download(sampled_fc$size()$getInfo()))
  dropped_coverage <- nrow(df_clean) - after_coverage
  if (dropped_coverage > 0) {
    sdm_warn(sprintf("%d point%s discarded (no satellite coverage, e.g. over water)",
                     dropped_coverage, if (dropped_coverage == 1) "" else "s"), indent = 1L)
  }

  # Deduplicate on embedding vectors — nearby points in the same pixel return
  # identical embeddings and add no new information to the classifier.
  emb_cols_r <- as.list(sprintf("A%02d", 0:63))
  sampled_fc  <- sampled_fc$distinct(emb_cols_r)
  final_count <- as.numeric(retry_curl_download(sampled_fc$size()$getInfo()))
  dropped_emb <- after_coverage - final_count
  if (dropped_emb > 0) {
    sdm_info(sprintf("%d point%s removed (duplicate embedding — same pixel or identical landscape)",
                     dropped_emb, if (dropped_emb == 1) "" else "s"), indent = 1L)
  }

  sdm_progress_done(pb_up)

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
    } else if (gee_method == "libsvm") {
      core_params <- c("numberOfTrees", "minLeafPopulation", "bagFraction", "shrinkage", "maxNodes", "variablesPerSplit", "lambda_", "polynomial", "batch_size")
      filtered_params <- params[setdiff(names(params), core_params)]
      # Default to linear kernel — RBF is O(n²) and times out on large training sets;
      # linear is sufficient for embedding spaces and scales with O(n×d).
      if (is.null(filtered_params$kernelType)) filtered_params$kernelType <- "LINEAR"
    } else {
      core_params <- c("numberOfTrees", "minLeafPopulation", "bagFraction", "shrinkage", "maxNodes", "variablesPerSplit", "lambda_", "polynomial", "batch_size")
      filtered_params <- params[setdiff(names(params), core_params)]
    }

    filtered_params <- filtered_params[!sapply(filtered_params, is.null)]
    clf <- do.call(clf_factory, filtered_params)
    clf <- clf$setOutputMode("PROBABILITY")

    trained_model <- clf$train(
      features = sampled_fc,
      classProperty = LABEL_COL,
      inputProperties = emb_cols
    )

    return(list(
      trained       = trained_model,
      is_classifier = TRUE,
      method        = method
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
      weights       = weights,
      intercept     = 0.0,
      is_classifier = FALSE,
      method        = method
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
