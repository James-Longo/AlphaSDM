#' Detect a Google Earth Engine synchronous-compute limit error
#'
#' Earth Engine's interactive `getInfo()` has a compute-time budget and a
#' 5000-feature page cap. Large point sets, many classes and expensive classifiers
#' exceed that budget and surface as one of the messages matched here. These are
#' the cases the batch-export path exists for.
#'
#' @param e A condition or a string.
#' @return TRUE when the message looks like a synchronous-compute limit.
#' @noRd
is_gee_timeout <- function(e) {
  msg <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
  grepl("Computation timed out", msg, fixed = TRUE) ||
    grepl("User memory limit exceeded", msg, fixed = TRUE) ||
    grepl("Collection query aborted", msg, fixed = TRUE) ||
    grepl("Too many concurrent aggregations", msg, fixed = TRUE) ||
    grepl("computation took too long", msg, ignore.case = TRUE)
}

#' Resolve a writable Earth Engine folder for temporary exports
#' @param project Earth Engine project id, or NULL to use the saved one.
#' @return An asset folder path. Errors when no writable folder can be found.
#' @noRd
async_asset_root <- function(project = NULL) {
  ee <- reticulate::import("ee")
  if (is.null(project) || !nzchar(project)) project <- .read_saved_project()
  if (!is.null(project) && nzchar(project)) return(sprintf("projects/%s/assets", project))

  # Last resort. On a Cloud project ee.data.getAssetRoots() returns the project's
  # assets rather than its roots, so the first entry is typically an image or table,
  # not a folder. Writing under it fails with "is neither a folder nor an image
  # collection", so only accept an entry that really is a folder.
  roots <- try(ee$data$getAssetRoots(), silent = TRUE)
  if (!inherits(roots, "try-error")) {
    for (r in roots) {
      if (identical(r[["type"]], "FOLDER") && !is.null(r[["id"]])) return(r[["id"]])
    }
  }
  stop("Async export needs a writable GEE asset folder. Pass gee_project = '<your-project>' ",
       "or run setup_gee() so the project id is saved.", call. = FALSE)
}

#' Read a FeatureCollection in pages, avoiding the 5000-feature cap
#'
#' Intended for a materialised asset. The values are already computed, so each page
#' is cheap and stays inside the synchronous compute limit. Paging a lazy collection
#' instead re-evaluates its graph once per page.
#'
#' @param fc An `ee$FeatureCollection`.
#' @param page_size Features per request.
#' @return A list with a `features` element, matching the shape `getInfo()` returns.
#' @noRd
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

#' Start an Export.table.toAsset batch task without waiting for it
#'
#' Start several exports before awaiting any of them and the chunks run
#' concurrently server-side.
#'
#' @param fc An `ee$FeatureCollection` to export.
#' @param project Earth Engine project id, or NULL to use the saved one.
#' @param select Optional property names to keep.
#' @return A list with the running `task` and its `asset_id`.
#' @noRd
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

#' Poll an export task until it finishes
#'
#' Prints roughly one line a minute so a long export shows progress rather than
#' sitting silent.
#'
#' @param handle A handle from `ee_start_fc_export()`.
#' @param poll_seconds Seconds between status checks.
#' @param max_minutes Overall limit. NULL, the default, means wait for as long as
#'   the task keeps running. A deadline on work that is progressing only turns
#'   Earth Engine's success into our own failure, which is what a fixed limit here
#'   kept doing. Earth Engine ends its own tasks, so this does not wait forever.
#'   Set ALPHASDM_MAX_WAIT_MINUTES for an unattended run that must not block.
#' @param max_queue_minutes Give up if the task has not started within this long.
#'   NULL, the default, waits. Queueing is how Earth Engine schedules work, not a
#'   sign that something is wrong, and a task that has not started yet has cost
#'   nothing to keep waiting for. Provided for an unattended run that would rather
#'   fall back than sit in a queue.
#' @return The asset id. Errors when the task fails, is cancelled, or gives up.
#' @noRd
ee_await_export <- function(handle, poll_seconds = 15, max_minutes = NULL,
                            max_queue_minutes = NULL) {
  env_cap <- suppressWarnings(as.numeric(Sys.getenv("ALPHASDM_MAX_WAIT_MINUTES", "")))
  if (!is.na(env_cap) && env_cap > 0) max_minutes <- env_cap
  deadline <- if (is.null(max_minutes)) NULL else Sys.time() + max_minutes * 60
  queue_deadline <- if (is.null(max_queue_minutes)) NULL else Sys.time() + max_queue_minutes * 60
  ever_ran <- FALSE
  start <- Sys.time(); last_beat <- 0; last_state <- ""
  repeat {
    st    <- handle$task$status()
    state <- st[["state"]]
    if (identical(state, "COMPLETED")) break
    if (state %in% c("FAILED", "CANCELLED", "CANCEL_REQUESTED")) {
      msg <- if (!is.null(st[["error_message"]])) st[["error_message"]] else state
      stop(sprintf("Async GEE export %s: %s", state, msg))
    }
    if (identical(state, "RUNNING")) ever_ran <- TRUE
    # Abandon a task that never got off the queue, but keep waiting on one that is
    # running: the wait that is worth cutting short is the backlog, not the work.
    if (!ever_ran && !is.null(queue_deadline) && Sys.time() > queue_deadline) {
      try(handle$task$cancel(), silent = TRUE)
      stop(sprintf("Async GEE export still queued after %g min (asset %s)",
                   max_queue_minutes, handle$asset_id))
    }
    if (!is.null(deadline) && Sys.time() > deadline) {
      try(handle$task$cancel(), silent = TRUE)
      stop(sprintf("Async GEE export exceeded max_minutes = %g (asset %s)",
                   max_minutes, handle$asset_id))
    }
    # Report the task's own state on every change, plus a heartbeat about once a
    # minute. READY means queued and RUNNING means working, and on a throttled tier
    # the difference between them is most of the wait.
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    if (!identical(state, last_state) || elapsed - last_beat >= 60) {
      lbl <- switch(state, READY = "queued", RUNNING = "running", tolower(state))
      sdm_info(sprintf("export %s server-side ... (%s elapsed)", lbl,
                       if (elapsed < 600) sprintf("%.0fs", elapsed)
                       else sprintf("%.0f min", elapsed / 60)), indent = 2L)
      last_beat <- elapsed; last_state <- state
    }
    Sys.sleep(poll_seconds)
  }
  handle$asset_id
}

#' Show the status of recent AlphaSDM Earth Engine batch tasks
#'
#' Lists the export and classifier tasks AlphaSDM has started, with each task's state
#' and age. Call it from a second session to see what Earth Engine is doing while a
#' long run is in progress.
#'
#' @param active_only TRUE shows only pending and running tasks. FALSE also lists
#'   recently finished ones.
#' @param since_minutes Only include tasks created within this many minutes.
#' @return A data frame of tasks with `description`, `state` and `age_min`,
#'   invisibly. Also prints them.
#' @export
sdm_gee_status <- function(active_only = TRUE, since_minutes = 180) {
  ee  <- reticulate::import("ee")
  ops <- tryCatch(ee$data$listOperations(), error = function(e) NULL)
  if (is.null(ops) || length(ops) == 0) { sdm_info("No Earth Engine tasks found."); return(invisible(NULL)) }
  df <- do.call(rbind, lapply(ops, function(o) {
    m <- o$metadata
    data.frame(description = tryCatch(m$description, error = function(e) NA_character_),
               state       = tryCatch(m$state,       error = function(e) NA_character_),
               created     = tryCatch(m$createTime,  error = function(e) NA_character_),
               stringsAsFactors = FALSE)
  }))
  df$created <- as.POSIXct(df$created, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
  df$age_min <- round(as.numeric(Sys.time() - df$created, units = "mins"), 1)
  keep <- !is.na(df$created) & df$age_min <= since_minutes & grepl("alphasdm", df$description)
  if (active_only) keep <- keep & df$state %in% c("PENDING", "RUNNING")
  df <- df[keep, , drop = FALSE]
  df <- df[order(df$created), , drop = FALSE]
  if (nrow(df) == 0) { sdm_info("No matching AlphaSDM Earth Engine tasks."); return(invisible(df)) }
  sdm_section(sprintf("AlphaSDM Earth Engine tasks (%d)", nrow(df)))
  for (i in seq_len(nrow(df)))
    sdm_info(sprintf("[%-9s] %s  (%.1f min ago)", df$state[i], df$description[i], df$age_min[i]), indent = 1L)
  invisible(df)
}

#' Export a FeatureCollection to a temporary Earth Engine asset and wait for it
#'
#' A batch task has no synchronous compute limit, so this path carries work that
#' `getInfo()` cannot.
#'
#' WARNING: this writes an asset and does not remove it. The caller has to delete
#' it with `ee_delete_asset_quietly()`.
#'
#' @param fc An `ee$FeatureCollection` to export.
#' @param project Earth Engine project id, or NULL to use the saved one.
#' @param select Optional property names to keep.
#' @param poll_seconds Seconds between status checks.
#' @param max_minutes Give up after this long.
#' @return The new asset id.
#' @noRd
ee_export_fc_to_asset <- function(fc, project = NULL, select = NULL,
                                  poll_seconds = 15, max_minutes = NULL) {
  handle <- ee_start_fc_export(fc, project, select)
  ee_await_export(handle, poll_seconds, max_minutes)
}

#' Delete an Earth Engine asset, ignoring failure
#' @param asset_id Asset to delete.
#' @return Nothing. Deletes the asset if it exists.
#' @noRd
ee_delete_asset_quietly <- function(asset_id) {
  ee <- reticulate::import("ee")
  try(ee$data$deleteAsset(asset_id), silent = TRUE)
  invisible(NULL)
}

#' Compute a FeatureCollection once into an asset and return it
#'
#' An Earth Engine collection is lazy, so a `sampleRegions` over the 64 embedding
#' bands is re-run on every request that touches it. Exporting it once to an asset
#' means training, classifying and scoring all read stored values instead. This is
#' what keeps large jobs inside the compute and memory limits.
#'
#' The data stays on Earth Engine throughout.
#'
#' WARNING: this writes an asset and does not remove it. The caller has to delete
#' it with `ee_delete_asset_quietly()`.
#'
#' @param fc An `ee$FeatureCollection` to materialise.
#' @param project Earth Engine project id, or NULL to use the saved one.
#' @param select Optional property names to keep.
#' @param poll_seconds Seconds between status checks.
#' @param max_minutes Give up after this long.
#' @return A list with the materialised `fc` and its `asset_id`.
#' @noRd
ee_materialize_fc_async <- function(fc, project = NULL, select = NULL,
                                    poll_seconds = 15, max_minutes = NULL) {
  ee <- reticulate::import("ee")
  asset_id <- ee_export_fc_to_asset(fc, project, select, poll_seconds, max_minutes)
  list(fc = ee$FeatureCollection(asset_id), asset_id = asset_id)
}

#' Store a trained classifier in an asset and load it back
#'
#' Used inline, a large tree model lives as one value in the same graph as
#' `classify`, and deep or many-class forests fail there with "Computed value is too
#' large". `Export.classifier.toAsset` trains as its own batch task and writes the
#' model to an asset, and `ee.Classifier.load()` then refers to it as stored data, so
#' the classify graph no longer carries the model itself.
#'
#' WARNING: this writes an asset and does not remove it. The caller has to delete
#' it with `ee_delete_asset_quietly()`.
#'
#' @param clf A trained `ee$Classifier`.
#' @param project Earth Engine project id, or NULL to use the saved one.
#' @param poll_seconds Seconds between status checks.
#' @param max_minutes Give up after this long.
#' @param max_queue_minutes Give up if the task has not started within this long.
#' @return A list with the reloaded `classifier` and its `asset_id`.
#' @noRd
ee_persist_classifier <- function(clf, project = NULL, poll_seconds = 15, max_minutes = NULL,
                                  max_queue_minutes = NULL) {
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
  ee_await_export(list(task = task, asset_id = asset_id), poll_seconds, max_minutes,
                  max_queue_minutes = max_queue_minutes)
  list(classifier = ee$Classifier$load(asset_id), asset_id = asset_id)
}

#' Remove leftover AlphaSDM temporary assets
#'
#' The package writes temporary Earth Engine assets while it works and deletes them
#' when it finishes. A run that is killed, crashes, or loses its connection never
#' reaches that cleanup, so the asset stays and counts against the project's storage
#' quota. This removes those leftovers.
#'
#' Only assets this package created are considered, matched on the `alphasdm_` name
#' it gives them. An asset is kept if a batch task for it is still pending or
#' running, or if it is newer than `older_than_hours`, so a job in progress in
#' another session is not disturbed. A large export can run for many hours, which is
#' why the default is deliberately generous.
#'
#' @param older_than_hours Keep assets younger than this. Default 48.
#' @param dry_run If TRUE, report what would be deleted and delete nothing.
#' @param project Earth Engine project id, or NULL to use the saved one.
#' @param quiet If TRUE, do not print anything.
#' @return The asset ids removed, or the ones that would be, invisibly.
#' @export
sdm_clean_assets <- function(older_than_hours = 48, dry_run = FALSE,
                             project = NULL, quiet = FALSE) {
  ee <- reticulate::import("ee")
  ensure_gee_authenticated(project = project)
  root <- async_asset_root(project)

  assets <- tryCatch(ee$data$listAssets(list(parent = root))[["assets"]],
                     error = function(e) NULL)
  if (is.null(assets) || length(assets) == 0) {
    if (!quiet) sdm_info("No Earth Engine assets found.")
    return(invisible(character(0)))
  }
  ids  <- vapply(assets, function(a) a[["id"]], character(1))
  mine <- grep("alphasdm_(async|clf|img)_[0-9]{14}_", ids, value = TRUE)
  if (length(mine) == 0) {
    if (!quiet) sdm_info("No leftover AlphaSDM assets.")
    return(invisible(character(0)))
  }

  # Never touch an asset whose task is still queued or running.
  active <- tryCatch({
    ops <- ee$data$listOperations()
    st  <- vapply(ops, function(o) tryCatch(o$metadata$state, error = function(e) ""), character(1))
    ds  <- vapply(ops, function(o) tryCatch(o$metadata$description, error = function(e) ""), character(1))
    ds[st %in% c("PENDING", "RUNNING")]
  }, error = function(e) character(0))

  stamp <- as.POSIXct(sub(".*_([0-9]{14})_.*", "\\1", mine), format = "%Y%m%d%H%M%S", tz = "")
  age_h <- as.numeric(difftime(Sys.time(), stamp, units = "hours"))
  suffix <- sub(".*_[0-9]{14}_", "", mine)

  drop <- !is.na(age_h) & age_h > older_than_hours &
    !vapply(suffix, function(sfx) any(grepl(sfx, active, fixed = TRUE)), logical(1))
  targets <- mine[drop]

  if (length(targets) == 0) {
    if (!quiet) sdm_info(sprintf("%d AlphaSDM asset(s) present, none old enough to remove.",
                                 length(mine)))
    return(invisible(character(0)))
  }
  if (!dry_run) for (a in targets) ee_delete_asset_quietly(a)
  if (!quiet) {
    sdm_info(sprintf("%s %d leftover AlphaSDM asset%s (older than %g h).",
                     if (dry_run) "Would remove" else "Removed", length(targets),
                     if (length(targets) == 1) "" else "s", older_than_hours))
  }
  invisible(targets)
}

#' Compute an image once into an asset and load it back
#'
#' The interactive endpoint that serves tile downloads has a small per-request
#' memory budget, and an image carrying a fitted model exceeds it at fine scales.
#' Earth Engine's own guidance for that case is to run the computation through the
#' batch system instead, which has a larger memory allowance and a longer timeout.
#' The result is a stored image, so reading tiles from it costs a read rather than
#' re-running the model over every pixel.
#'
#' Progress is reported to the console as the task runs, so the caller does not need
#' to watch the Earth Engine task list.
#'
#' WARNING: this writes an asset and does not remove it. The caller has to delete it
#' with `ee_delete_asset_quietly()`.
#'
#' @param img An `ee$Image` to compute.
#' @param region Geometry to compute over.
#' @param scale Pixel size in metres.
#' @param project Earth Engine project id, or NULL to use the saved one.
#' @param poll_seconds Seconds between status checks.
#' @param max_minutes Give up after this long.
#' @return A list with the stored `img` and its `asset_id`.
#' @noRd
ee_persist_image <- function(img, region, scale, project = NULL,
                             poll_seconds = 15, max_minutes = NULL) {
  ee     <- reticulate::import("ee")
  root   <- async_asset_root(project)
  stamp  <- format(Sys.time(), "%Y%m%d%H%M%S")
  suffix <- paste(sample(c(letters, 0:9), 6, replace = TRUE), collapse = "")
  asset_id <- sprintf("%s/alphasdm_img_%s_%s", root, stamp, suffix)

  task <- ee$batch$Export$image$toAsset(
    image = img, description = paste0("alphasdm_img_", suffix), assetId = asset_id,
    region = region, scale = scale, maxPixels = 1e13
  )
  task$start()
  sdm_info(sprintf("Computing the map server-side -> %s", basename(asset_id)), indent = 2L)
  ee_await_export(list(task = task, asset_id = asset_id), poll_seconds, max_minutes)
  list(img = ee$Image(asset_id), asset_id = asset_id)
}

#' Read a FeatureCollection through a batch export instead of getInfo()
#'
#' Export to a temporary asset, read it back in pages, then delete it. Slower than
#' `getInfo()` but under no synchronous compute limit, so it carries work that
#' `getInfo()` refuses.
#'
#' Unlike the other export helpers here, this one deletes the asset it creates.
#'
#' @param fc An `ee$FeatureCollection` to read.
#' @param project Earth Engine project id, or NULL to use the saved one.
#' @param select Optional property names to keep.
#' @param poll_seconds Seconds between status checks.
#' @param max_minutes Give up after this long.
#' @return A list with a `features` element, matching the shape `getInfo()` returns.
#' @noRd
ee_table_to_info_async <- function(fc, project = NULL, select = NULL,
                                   poll_seconds = 15, max_minutes = NULL) {
  ee <- reticulate::import("ee")
  asset_id <- ee_export_fc_to_asset(fc, project, select, poll_seconds, max_minutes)
  info <- read_fc_paged(ee$FeatureCollection(asset_id))
  ee_delete_asset_quietly(asset_id)
  info
}
