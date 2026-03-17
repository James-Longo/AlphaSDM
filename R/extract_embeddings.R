#' Internal helper to retrieve properties from a FeatureCollection in chunks
#' @keywords internal
retrieve_gee_properties_chunked <- function(fc, selectors, chunk_size = 1000) {
    ee <- reticulate::import("ee")
    label <- "Retrieving GEE results"
    t0 <- log_time(label)
    
    n_available <- as.numeric(fc$size()$getInfo())
    if (n_available == 0) return(list())
    
    full_list <- fc$toList(n_available)
    res_list <- list()
    for (i in seq(0, n_available - 1, by = chunk_size)) {
        start_idx <- i
        end_idx <- min(i + chunk_size, n_available)
        
        # Professional progress logging
        if (n_available > chunk_size) {
            message(sprintf("  Retrieving chunk %i/%i...", 
                            (i %/% chunk_size) + 1, 
                            ceiling(n_available / chunk_size)))
        }
        
        chunk_fc <- ee$FeatureCollection(full_list$slice(start_idx, end_idx))
        
        chunk_res <- chunk_fc$reduceColumns(
            reducer = ee$Reducer$toList()$`repeat`(length(selectors)),
            selectors = as.list(selectors)
        )$getInfo()
        
        if (!is.null(chunk_res$list)) {
            res_list <- c(res_list, chunk_res$list)
        }
    }
    
    log_time(label, t0)
    return(res_list)
}

#' Extract Alpha Earth Embeddings
#'
#' This function extracts the 64-dimensional Alpha Earth satellite embeddings for a set of coordinates.
#' Standardizes extraction using Google Earth Engine.
#'
#' @param df A data frame containing at least `longitude`, `latitude`, and `year`.
#' @param scale Optional. Resolution in meters for extraction. Defaults to 10.
#' @param gee_project Optional. Google Cloud Project ID. If NULL, uses local config or environment variables.
#' @return A data frame containing the original columns plus 64 embedding columns (A00-A63).
#' @export
extract_embeddings <- function(df, scale = 10, gee_project = NULL) {
    # 1. GEE Auth
    ensure_gee_authenticated(project = gee_project)

    # 2. Preparation and Sampling
    t_prep <- log_time(sprintf("Preparing %d points for GEE extraction", nrow(df)))
    
    ee <- reticulate::import("ee")
    
    # Create individual points to retain order/ID
    # Chunked to be safe with memory and R-Python bridge
    chunk_size <- 2000
    all_features <- list()
    for (i in seq(1, nrow(df), by = chunk_size)) {
        end <- min(i + chunk_size - 1, nrow(df))
        chunk_df <- df[i:end, ]
        chunk_features <- lapply(seq_len(nrow(chunk_df)), function(j) {
            ee$Feature(
                ee$Geometry$Point(c(as.numeric(chunk_df$longitude[j]), as.numeric(chunk_df$latitude[j]))),
                list(year = as.integer(chunk_df$year[j]), tmp_id = as.integer(i + j - 2))
            )
        })
        all_features <- c(all_features, chunk_features)
    }
    upload_fc <- ee$FeatureCollection(all_features)
    log_time("Point preparation", t_prep)
    
    # Sample embeddings
    t_samp <- log_time("Sampling Alpha Earth embeddings on GEE")
    sampled_fc <- get_embeddings_at_fc(upload_fc, scale, properties = c("year", "tmp_id"))
    log_time("Sampling", t_samp)
    
    # 3. Retrieve results to local data frame
    emb_cols <- sprintf("A%02d", 0:63)
    selectors <- c(emb_cols, "year", "tmp_id")
    
    res_list <- retrieve_gee_properties_chunked(sampled_fc, selectors)
    
    # 4. Align with input data (Scientific Integrity)
    t_df <- log_time("Aligning results with input data")
    
    # Create template with all original IDs
    final_df <- data.frame(tmp_id = 0:(nrow(df) - 1))
    
    if (length(res_list) > 0) {
        res_df <- as.data.frame(do.call(rbind, res_list))
        colnames(res_df) <- selectors
        
        # Ensure numeric for merging
        res_df$tmp_id <- as.numeric(res_df$tmp_id)
        for (col in emb_cols) res_df[[col]] <- as.numeric(res_df[[col]])
        
        # Left join to preserve all original indices
        final_df <- merge(final_df, res_df, by = "tmp_id", all.x = TRUE)
    } else {
        # Fill with NAs if nothing returned
        for (col in selectors[selectors != "tmp_id"]) final_df[[col]] <- NA
    }
    
    # Reiterate original spatial info for absolute traceability
    final_df$longitude <- df$longitude
    final_df$latitude <- df$latitude
    
    # Restore original order (merge might sort)
    final_df <- final_df[order(final_df$tmp_id), ]
    final_df$tmp_id <- NULL
    
    # Log scientific summary
    n_missing <- sum(is.na(final_df$A00))
    if (n_missing > 0) {
        message(sprintf("[autoSDM] Scientific Alert: %d points (%.1f%%) were outside satellite coverage/masked.", 
                        n_missing, (n_missing / nrow(df)) * 100))
    }
    
    log_time("Data alignment", t_df)
    return(final_df)
}

#' Initialize Earth Engine with Service Account
#'
#' @param json_path Path to the service account JSON key. Defaults to GEE_SERVICE_ACCOUNT_KEY env var.
#' @param venv_path Path to the Python virtual environment.
#' @export
#' @keywords internal
ee_auth_service <- function(json_path = Sys.getenv("GEE_SERVICE_ACCOUNT_KEY"), venv_path = NULL) {
    if (json_path == "" || !file.exists(json_path)) {
        stop("Valid service account JSON key path required.")
    }

    # We can't easily switch projects once ee is initialized in reticulate without restarting,
    # but we can try to rely on the CLI's --key argument for subsequent calls.
    # This function primarily exists to register the key path for the session.
    assign("sa_json_key", json_path, envir = .GlobalEnv)
    message("Service account registered for this session.")
}
