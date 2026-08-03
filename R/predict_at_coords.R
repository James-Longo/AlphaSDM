#' Score coordinates with fitted AlphaSDM models
#'
#' Samples the Alpha Earth embeddings at a set of coordinates and scores them with
#' one or more already-fitted models, entirely server-side. Results are pulled back
#' in parallel CSV batches. Use [evaluate_models()] when you also want cross-validation
#' or accuracy metrics; this function is the lighter path for scoring new points with
#' models you already have in hand.
#'
#' @param df Data frame with columns `longitude`, `latitude` and `year`.
#' @param models Named list of fitted model objects, as returned in the `models`
#'   element of the internal training step (each carries `method`, `is_classifier`
#'   and either a trained GEE classifier or centroid `weights`). Names are used to
#'   label the output columns.
#' @param scale Embedding resolution in metres.
#' @param aoi_year Year of the Alpha Earth mosaic to sample.
#' @param gee_project Optional Earth Engine project override (normally set via [setup_gee()]).
#' @param options Named list of advanced options; `batch_size` sets the number of
#'   points per download batch.
#' @return A list with `data` (the input data frame plus a `similarity` column, or
#'   `similarity_<model>` columns when several models are supplied) and `metrics`.
#' @export
predict_at_coords <- function(df, models, scale = NULL, aoi_year = NULL,
                              gee_project = NULL, options = list()) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)

  if (is.null(scale)) stop("Parameter 'scale' must be provided.")
  if (is.null(aoi_year)) stop("Parameter 'aoi_year' must be provided.")

  sdm_section("Predicting at coordinates (server-side)")

  # Row key carried through GEE so scores can be matched back to the input rows
  # regardless of the order features come back in.
  df$tmp_id <- seq_len(nrow(df)) - 1L

  # Upload points and sample embeddings (both stay on GEE).
  upload_fc  <- upload_points_to_gee(df)
  sampled_fc <- get_embeddings_at_fc(upload_fc, scale, properties = c("year", "tmp_id"))

  final_results <- df

  sdm_info(sprintf("Scoring with %d model%s",
                   length(models), if (length(models) == 1) "" else "s"))

  models_to_score <- list()
  for (i in seq_along(models)) {
    model_res <- models[[i]]
    if (!is.list(model_res)) next
    m_name <- names(models)[i]
    if (is.null(m_name) || !nzchar(m_name)) m_name <- model_res$method
    models_to_score[[m_name]] <- model_res
  }
  if (length(models_to_score) == 0) stop("No fitted models supplied to score with.")

  scored_fc <- predict_all_models_gee(sampled_fc, models_to_score)

  # Download results in parallel batches.
  selectors <- as.list(c("tmp_id", paste0("pred_", names(models_to_score))))

  chunk_size <- if (!is.null(options$batch_size)) {
    as.integer(options$batch_size)
  } else {
    safe_chunk_size <- if (scale <= 160) 5000 else 10000
    max(safe_chunk_size, ceiling(nrow(df) / 40))
  }

  batch_starts <- seq(1, nrow(df), by = chunk_size)

  sdm_info(sprintf("Downloading results across %d parallel batch%s ...",
                   length(batch_starts),
                   if (length(batch_starts) == 1) "" else "es"))

  process_batch <- function(idx) {
    # Re-import ee inside the parallel worker to avoid closure/fork issues.
    ee      <- reticulate::import("ee")
    i       <- batch_starts[idx]
    end_idx <- min(i + chunk_size - 1, nrow(df))
    chunk_fc <- scored_fc$filter(ee$Filter$gte("tmp_id", as.integer(i - 1)))$
                          filter(ee$Filter$lte("tmp_id", as.integer(end_idx - 1)))

    for (attempt in seq_len(5L)) {
      res <- tryCatch({
        url <- chunk_fc$getDownloadURL(filetype = "CSV", selectors = selectors)
        read.csv(url)
      }, error = function(e) e)
      if (is.data.frame(res)) return(res)
      Sys.sleep(2^(attempt - 1) + runif(1, 0, 1))
    }
    data.frame()
  }

  all_res_list <- future_lapply(seq_along(batch_starts), process_batch, future.seed = TRUE)
  res_df       <- do.call(rbind, all_res_list)

  if (!is.null(res_df) && nrow(res_df) > 0) {
    res_df     <- res_df[order(res_df$tmp_id), ]
    idx        <- match(res_df$tmp_id, (seq_len(nrow(df)) - 1))
    valid_mask <- !is.na(idx)
    for (m in names(models_to_score)) {
      col_name <- if (length(models_to_score) > 1) paste0("similarity_", m) else "similarity"
      final_results[idx[valid_mask], col_name] <- as.numeric(res_df[[paste0("pred_", m)]][valid_mask])
    }
  }

  final_results$tmp_id <- NULL
  list(data = final_results, metrics = list())
}
