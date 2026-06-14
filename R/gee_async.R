#' Detect a Google Earth Engine synchronous-compute limit error
#'
#' GEE's interactive `getInfo()` has a compute-time budget (and a 5000-feature
#' page cap). Heavy jobs — large point sets, many classes, expensive classifiers —
#' blow that budget and surface as one of the messages matched here. These are
#' exactly the cases the async (batch-export) path is designed for.
#' @keywords internal
is_gee_timeout <- function(e) {
  msg <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
  grepl("Computation timed out", msg, fixed = TRUE) ||
    grepl("User memory limit exceeded", msg, fixed = TRUE) ||
    grepl("Collection query aborted", msg, fixed = TRUE) ||
    grepl("Too many concurrent aggregations", msg, fixed = TRUE) ||
    grepl("computation took too long", msg, ignore.case = TRUE)
}

#' Standard message recommending the async path after a synchronous timeout
#' @keywords internal
async_timeout_message <- function(e, where = "scoring") {
  paste0(
    "GEE timed out during ", where, " (", conditionMessage(e), ").\n",
    "  This job is too large/complex for synchronous evaluation. Re-run with ",
    "async = TRUE\n  to materialise results via a GEE batch export (no synchronous ",
    "compute limit)."
  )
}

#' Resolve a writable GEE asset folder for temporary async exports
#' @keywords internal
async_asset_root <- function(project = NULL) {
  ee <- reticulate::import("ee")
  if (!is.null(project) && nzchar(project)) return(sprintf("projects/%s/assets", project))
  roots <- try(ee$data$getAssetRoots(), silent = TRUE)
  if (!inherits(roots, "try-error") && length(roots) > 0 && !is.null(roots[[1]][["id"]])) {
    return(roots[[1]][["id"]])
  }
  stop("Async export needs a writable GEE asset folder. Pass gee_project = '<your-project>'.")
}

#' Read a FeatureCollection in pages via computeFeatures (no 5000-feature cap)
#'
#' Used to read a materialised asset; because the asset is already computed, each
#' page is cheap and won't hit the synchronous compute limit.
#' @keywords internal
read_fc_paged <- function(fc, page_size = 5000L) {
  ee_data <- reticulate::import("ee.data")
  feats <- list(); tok <- NULL
  repeat {
    params <- list(expression = fc, pageSize = as.integer(page_size))
    if (!is.null(tok) && nzchar(tok)) params$pageToken <- tok
    r <- retry_curl_download(ee_data$computeFeatures(params))
    if (!is.null(r[["features"]])) feats <- c(feats, r[["features"]])
    tok <- r[["nextPageToken"]]
    if (is.null(tok) || !nzchar(tok)) break
  }
  list(features = feats)
}

#' Start (without waiting for) an Export.table.toAsset batch task
#'
#' Returns the running task handle and its target asset id. Starting many exports
#' before awaiting any of them lets independent chunks run concurrently server-side.
#' @keywords internal
ee_start_fc_export <- function(fc, project = NULL, select = NULL) {
  ee <- reticulate::import("ee")
  if (!is.null(select)) fc <- fc$select(as.list(select))

  root   <- async_asset_root(project)
  stamp  <- format(Sys.time(), "%Y%m%d%H%M%S")
  suffix <- paste(sample(c(letters, 0:9), 6, replace = TRUE), collapse = "")
  asset_id <- sprintf("%s/alphasdm_async_%s_%s", root, stamp, suffix)

  task <- ee$batch$Export$table$toAsset(
    collection = fc, description = paste0("alphasdm_async_", suffix), assetId = asset_id
  )
  task$start()
  sdm_info(sprintf("Async batch export started (server-side) -> %s", asset_id), indent = 1L)
  list(task = task, asset_id = asset_id)
}

#' Poll a started export task until it completes (or fail/timeout)
#' @keywords internal
ee_await_export <- function(handle, poll_seconds = 15, max_minutes = 60) {
  deadline <- Sys.time() + max_minutes * 60
  repeat {
    st    <- handle$task$status()
    state <- st[["state"]]
    if (identical(state, "COMPLETED")) break
    if (state %in% c("FAILED", "CANCELLED", "CANCEL_REQUESTED")) {
      msg <- if (!is.null(st[["error_message"]])) st[["error_message"]] else state
      stop(sprintf("Async GEE export %s: %s", state, msg))
    }
    if (Sys.time() > deadline) {
      try(handle$task$cancel(), silent = TRUE)
      stop(sprintf("Async GEE export exceeded max_minutes = %d (asset %s)", max_minutes, handle$asset_id))
    }
    Sys.sleep(poll_seconds)
  }
  handle$asset_id
}

#' Export a FeatureCollection to a temporary GEE asset (batch) and wait
#'
#' Runs an `Export.table.toAsset` batch task — server-side, with no synchronous
#' compute limit — and polls until it completes. Returns the new asset id. The
#' caller is responsible for deleting the asset (see `ee_delete_asset_quietly`).
#' @keywords internal
ee_export_fc_to_asset <- function(fc, project = NULL, select = NULL,
                                  poll_seconds = 15, max_minutes = 60) {
  handle <- ee_start_fc_export(fc, project, select)
  ee_await_export(handle, poll_seconds, max_minutes)
}

#' Delete a GEE asset, ignoring errors
#' @keywords internal
ee_delete_asset_quietly <- function(asset_id) {
  ee <- reticulate::import("ee")
  try(ee$data$deleteAsset(asset_id), silent = TRUE)
  invisible(NULL)
}

#' Materialise a FeatureCollection to a GEE asset and return it (kept server-side)
#'
#' The key primitive for heavy pipelines: a lazy FeatureCollection (e.g. a
#' `sampleRegions` of the 64-band embeddings at many points) is exported once to an
#' asset, so every downstream op (train, classify, score) reads PRE-COMPUTED values
#' instead of re-running the sampling/training graph per request. That is what
#' prevents the "Computation timed out" / "out of memory" failures on large,
#' many-class jobs. Returns the materialised `ee$FeatureCollection` plus its
#' `asset_id`; the caller should `ee_delete_asset_quietly(asset_id)` when done.
#' @keywords internal
ee_materialize_fc_async <- function(fc, project = NULL, select = NULL,
                                    poll_seconds = 15, max_minutes = 60) {
  ee <- reticulate::import("ee")
  asset_id <- ee_export_fc_to_asset(fc, project, select, poll_seconds, max_minutes)
  list(fc = ee$FeatureCollection(asset_id), asset_id = asset_id)
}

#' Classifier types that `Export.classifier.toAsset` supports
#' @keywords internal
PERSISTABLE_CLASSIFIERS <- c("rf", "cart")   # smileRandomForest, smileCart (+ decision trees)

#' Persist a trained classifier to a GEE asset and reload it
#'
#' Large tree models (deep / many-class / huge training sets) fail with "Computed
#' value is too large" when used inline, because the whole model lives as one value
#' in the same graph as `classify`. `Export.classifier.toAsset` runs the training as
#' its own batch task and writes the model to an asset; `ee.Classifier.load()` then
#' references it as stored data, so the classify graph no longer carries the oversized
#' training value. Returns the reloaded `ee$Classifier` and the `asset_id` to clean up
#' with `ee_delete_asset_quietly()`.
#' @keywords internal
ee_persist_classifier <- function(clf, project = NULL, poll_seconds = 15, max_minutes = 60) {
  ee     <- reticulate::import("ee")
  root   <- async_asset_root(project)
  stamp  <- format(Sys.time(), "%Y%m%d%H%M%S")
  suffix <- paste(sample(c(letters, 0:9), 6, replace = TRUE), collapse = "")
  asset_id <- sprintf("%s/alphasdm_clf_%s_%s", root, stamp, suffix)

  task <- ee$batch$Export$classifier$toAsset(
    classifier = clf, description = paste0("alphasdm_clf_", suffix), assetId = asset_id
  )
  task$start()
  sdm_info(sprintf("Persisting classifier via batch export -> %s", asset_id), indent = 1L)
  ee_await_export(list(task = task, asset_id = asset_id), poll_seconds, max_minutes)
  list(classifier = ee$Classifier$load(asset_id), asset_id = asset_id)
}

#' Sample 64-band embeddings at many points and materialise them to a GEE asset (chunked)
#'
#' Sampling + exporting a very large or continent-wide point set in one task blows
#' GEE's memory. This samples and exports in independent chunks (each light enough to
#' fit), starts the chunk exports concurrently, awaits them, then merges the per-chunk
#' assets into ONE materialised FeatureCollection. Downstream training/classification
#' then read pre-computed embeddings (no re-sampling), which is what makes heavy,
#' many-class jobs tractable. Embeddings never leave GEE. Returns the merged
#' `ee$FeatureCollection` and the `asset_ids` to clean up with `ee_delete_asset_quietly`.
#'
#' @param df          data frame with longitude, latitude, year (+ any `properties`).
#' @param scale       sampling scale (m).
#' @param properties  point properties to carry through (besides the A00–A63 bands).
#' @param years       list of years to sample (passed to `get_embeddings_at_fc`).
#' @param project     GEE project id for the temp asset folder.
#' @param chunk_size  points per chunk/export task.
#' @keywords internal
sample_embeddings_to_asset <- function(df, scale, properties, years, project = NULL,
                                       chunk_size = 4000L, poll_seconds = 15, max_minutes = 60) {
  ee  <- reticulate::import("ee")
  emb <- sprintf("A%02d", 0:63)
  starts <- seq(1L, nrow(df), by = as.integer(chunk_size))

  # Start every chunk's sample+export, then await them (concurrent server-side).
  handles <- lapply(starts, function(s) {
    sub    <- df[s:min(s + chunk_size - 1L, nrow(df)), , drop = FALSE]
    fc_sub <- get_embeddings_at_fc(upload_points_to_gee(sub), scale,
                                   properties = properties, geometries = TRUE, years = years)
    ee_start_fc_export(fc_sub$select(as.list(c(emb, properties))), project = project)
  })
  asset_ids <- vapply(handles, function(h) ee_await_export(h, poll_seconds, max_minutes), character(1))

  merged <- Reduce(function(a, b) a$merge(b),
                   lapply(asset_ids, function(a) ee$FeatureCollection(a)))
  list(fc = merged, asset_ids = asset_ids)
}

#' Materialise a FeatureCollection via a GEE batch export, then read it
#'
#' The async alternative to `fc$getInfo()`: export to a temporary asset (batch, no
#' synchronous compute limit), read the materialised asset back in pages (cheap), and
#' delete it. Returns a list with a `features` element, matching `getInfo()`'s shape.
#' @keywords internal
ee_table_to_info_async <- function(fc, project = NULL, select = NULL,
                                   poll_seconds = 15, max_minutes = 60) {
  ee <- reticulate::import("ee")
  asset_id <- ee_export_fc_to_asset(fc, project, select, poll_seconds, max_minutes)
  info <- read_fc_paged(ee$FeatureCollection(asset_id))
  ee_delete_asset_quietly(asset_id)
  info
}
