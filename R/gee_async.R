#' Detect a Google Earth Engine compute-limit error
#'
#' Earth Engine's interactive requests have a compute-time and memory budget.
#' Large point sets and expensive classifiers exceed it and fail with one of the
#' messages matched here. These are the cases the batch-export path exists for.
#'
#' @param e A condition or a string.
#' @return TRUE when the message looks like a compute limit.
#' @noRd
is_gee_timeout <- function(e) {
  msg <- if (inherits(e, "condition")) conditionMessage(e) else as.character(e)
  grepl(GEE_LIMIT_PATTERN, msg, ignore.case = TRUE)
}

#' Messages Earth Engine uses for a request that exceeded its compute budget
#' @noRd
GEE_LIMIT_PATTERN <- paste(c(
  "Computation timed out", "User memory limit exceeded", "Collection query aborted",
  "Too many concurrent aggregations", "computation took too long", "out of memory"),
  collapse = "|")

#' A fresh name for a temporary asset in the user's project
#'
#' Every temporary asset AlphaSDM writes is named `alphasdm_<kind>_<timestamp>_<code>`,
#' which is what gee_clean_assets() matches on.
#' @param kind Short label, such as "table" or "clf".
#' @return list(id = full asset id, description = task description).
#' @noRd
temp_asset_id <- function(kind, project = NULL) {
  project <- .resolve_project(project)
  if (is.null(project))
    stop("Batch exports need an Earth Engine project. Run setup_gee() or pass ",
         "gee_project.", call. = FALSE)
  code <- paste(sample(c(letters, 0:9), 6, replace = TRUE), collapse = "")
  name <- sprintf("alphasdm_%s_%s_%s", kind, format(Sys.time(), "%Y%m%d%H%M%S"), code)
  list(id = sprintf("projects/%s/assets/%s", project, name), description = name)
}

#' Read a FeatureCollection in pages, avoiding the 5000-feature cap
#'
#' Intended for a stored asset. The values are already computed, so each page is
#' cheap and stays inside the interactive compute limit. Paging a lazy collection
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

#' Describe how far along an Earth Engine batch task is
#'
#' The Operations API carries a progress fraction and compute consumed, which
#' `task$status()` does not. Either may be absent on a young task.
#'
#' @param op_name Operation name from `task$status()[["name"]]`.
#' @return A string to append to a progress line, empty when nothing is reported yet.
#' @noRd
ee_task_progress <- function(op_name) {
  if (is.null(op_name) || !nzchar(op_name)) return("")
  ee <- reticulate::import("ee")
  m <- tryCatch(ee$data$getOperation(op_name)$metadata, error = function(e) NULL)
  if (is.null(m)) return("")
  one <- function(x) {
    v <- suppressWarnings(as.numeric(x))
    if (length(v) != 1L || is.na(v)) NULL else v
  }
  bits <- character(0)
  # One decimal: a long export can spend an hour inside a single percent.
  pct <- one(m$progress)
  if (!is.null(pct)) bits <- c(bits, sprintf("%.1f%%", 100 * pct))
  eecu <- one(m$batchEecuUsageSeconds)
  if (!is.null(eecu)) bits <- c(bits, sprintf("%.0f EECU-s", eecu))
  if (!length(bits)) "" else paste0(" [", paste(bits, collapse = ", "), "]")
}

#' Print the live task-monitor links, once per session
#'
#' Batch tasks can sit in Google's queue for hours when the monthly compute
#' quota is exhausted; the wait is visible, and tasks can be cancelled, in the
#' Code Editor Tasks tab and the Cloud Console Earth Engine page.
#' @noRd
sdm_task_monitor_hint <- function(project = NULL) {
  if (isTRUE(.alphasdm_env$task_hint_printed)) return(invisible())
  .alphasdm_env$task_hint_printed <- TRUE
  project <- .resolve_project(project)
  sdm_info("Watch batch tasks live: https://code.earthengine.google.com/tasks",
           indent = 1L)
  if (!is.null(project))
    sdm_info(sprintf(
      "  or in Cloud Console: https://console.cloud.google.com/earth-engine/tasks?project=%s",
      project), indent = 1L)
  invisible()
}

#' Wait for Earth Engine batch tasks to finish
#'
#' Waits for as long as the tasks run: Earth Engine ends its own tasks, and
#' queueing is normal scheduling rather than a fault. Prints a line when the
#' state changes and about once a minute otherwise.
#'
#' @param tasks Named list of started `ee.batch.Task` objects.
#' @param poll_seconds Seconds between status checks.
#' @return Named character vector of error messages for the tasks that failed
#'   or were cancelled; empty when all completed. The caller decides what a
#'   failure means.
#' @noRd
ee_await_tasks <- function(tasks, poll_seconds = 15) {
  sdm_task_monitor_hint()
  state <- setNames(rep("READY", length(tasks)), names(tasks))
  error <- setNames(character(0), character(0))
  start <- Sys.time(); last_beat <- -Inf; last_line <- ""
  repeat {
    for (nm in names(tasks)[!state %in% c("COMPLETED", "FAILED", "CANCELLED")]) {
      st <- tasks[[nm]]$status()
      state[[nm]] <- if (st[["state"]] == "CANCEL_REQUESTED") "CANCELLED" else st[["state"]]
      if (state[[nm]] %in% c("FAILED", "CANCELLED"))
        error[[nm]] <- if (!is.null(st[["error_message"]])) st[["error_message"]] else state[[nm]]
    }
    if (all(state %in% c("COMPLETED", "FAILED", "CANCELLED"))) break
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    running <- names(tasks)[state == "RUNNING"]
    line <- if (length(tasks) == 1L) {
      switch(state[[1]], READY = "queued", RUNNING = "running", tolower(state[[1]]))
    } else {
      sprintf("%d of %d done, %d running, %d queued", sum(state == "COMPLETED"),
              length(tasks), length(running), sum(state == "READY"))
    }
    if (!identical(line, last_line) || elapsed - last_beat >= 60) {
      detail <- if (length(running))
        ee_task_progress(tasks[[running[1]]]$status()[["name"]]) else ""
      sdm_info(sprintf("export %s%s (%s elapsed)", line, detail,
                       if (elapsed < 600) sprintf("%.0fs", elapsed)
                       else sprintf("%.0f min", elapsed / 60)), indent = 2L)
      last_beat <- elapsed; last_line <- line
    }
    Sys.sleep(poll_seconds)
  }
  error
}

#' Compute FeatureCollections into temporary assets and read them as one
#'
#' A collection is lazy, so a `sampleRegions` over the 64 embedding bands is
#' re-run by every request that touches it. Exporting it once means training,
#' classifying and scoring all read stored values instead. Each collection is
#' exported as its own task, so the tasks run concurrently and each export's
#' graph stays small.
#'
#' WARNING: the assets are not removed here; the caller deletes them.
#'
#' @param fcs List of `ee$FeatureCollection`s. Exported features need geometry.
#' @return list(fc = the stored collections merged, asset_ids).
#' @noRd
ee_store_tables <- function(fcs, project = NULL) {
  ee <- reticulate::import("ee")
  ids <- lapply(seq_along(fcs), function(i) temp_asset_id("table", project))
  tasks <- lapply(seq_along(fcs), function(i) {
    task <- ee$batch$Export$table$toAsset(collection = fcs[[i]],
                                          description = ids[[i]]$description,
                                          assetId = ids[[i]]$id)
    task$start()
    task
  })
  names(tasks) <- vapply(ids, `[[`, "", "description")
  failed <- ee_await_tasks(tasks)
  asset_ids <- vapply(ids, `[[`, "", "id")
  if (length(failed)) {
    for (a in asset_ids) ee_delete_asset_quietly(a)
    stop(sprintf("Storing the table on Earth Engine failed: %s",
                 paste(unique(failed), collapse = "; ")), call. = FALSE)
  }
  stored <- lapply(asset_ids, ee$FeatureCollection)
  list(fc = Reduce(function(a, b) a$merge(b), stored), asset_ids = asset_ids)
}

#' Read a FeatureCollection through a batch export instead of getInfo()
#'
#' Slower than an interactive read but under no interactive compute limit, so it
#' carries work that `getInfo()` refuses. The temporary asset is deleted.
#'
#' @return A list with a `features` element, matching the shape `getInfo()` returns.
#' @noRd
ee_read_via_batch <- function(fc, project = NULL) {
  ee <- reticulate::import("ee")
  stored <- ee_store_tables(list(fc), project)
  on.exit(ee_delete_asset_quietly(stored$asset_ids), add = TRUE)
  read_fc_paged(stored$fc)
}

#' Store a trained classifier in an asset and load it back
#'
#' Inline, a large tree model is one value in the same graph as `classify`, and
#' big forests fail there with "Computed value is too large".
#' `Export.classifier.toAsset` writes the model to an asset, and
#' `ee.Classifier.load()` then refers to it as stored data. Earth Engine can
#' store only random forests and CART.
#'
#' WARNING: the asset is not removed here; the caller deletes it.
#'
#' @param clf A trained `ee$Classifier`.
#' @return A list with the reloaded `classifier` and its `asset_id`.
#' @noRd
ee_persist_classifier <- function(clf, project = NULL) {
  ee <- reticulate::import("ee")
  id <- temp_asset_id("clf", project)
  task <- ee$batch$Export$classifier$toAsset(classifier = clf, description = id$description,
                                             assetId = id$id)
  task$start()
  sdm_info(sprintf("Storing the classifier -> %s", id$id), indent = 1L)
  failed <- ee_await_tasks(setNames(list(task), id$description))
  if (length(failed)) stop(sprintf("Storing the classifier failed: %s", failed), call. = FALSE)
  list(classifier = ee$Classifier$load(id$id), asset_id = id$id)
}

#' Delete an Earth Engine asset, ignoring failure
#' @param asset_id Asset id(s) to delete.
#' @noRd
ee_delete_asset_quietly <- function(asset_id) {
  ee <- reticulate::import("ee")
  for (a in asset_id) try(ee$data$deleteAsset(a), silent = TRUE)
  invisible(NULL)
}

#' Show recent AlphaSDM Earth Engine tasks
#'
#' Lists the export tasks AlphaSDM has started, with each task's state and age.
#' Call it from a second R session to see what Earth Engine is doing while a
#' long run is in progress. To check the connection itself, use
#' [gee_status()].
#'
#' @param active_only TRUE shows only pending and running tasks; FALSE also
#'   lists recently finished ones.
#' @param since_minutes Only include tasks created within this many minutes.
#' @return A data frame of tasks with `description`, `state` and `age_min`,
#'   invisibly. Also prints them.
#' @examples
#' \dontrun{
#' gee_tasks(active_only = FALSE)
#' }
#' @export
gee_tasks <- function(active_only = TRUE, since_minutes = 180) {
  ensure_gee_authenticated()
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

#' Remove leftover AlphaSDM temporary assets
#'
#' AlphaSDM writes temporary Earth Engine assets while it works and deletes them
#' when it finishes. A run that is killed or loses its connection never reaches
#' that cleanup, so the asset stays and counts against the project's storage
#' quota. This removes those leftovers.
#'
#' Only assets named by AlphaSDM (`alphasdm_...`) are considered. An asset is
#' kept if a task for it is still pending or running, or if it is newer than
#' `older_than_hours`, so a job in progress in another session is not disturbed.
#'
#' @param older_than_hours Keep assets younger than this. Default 48.
#' @param dry_run If TRUE, report what would be deleted and delete nothing.
#' @param project Earth Engine project id, or NULL to use the saved one.
#' @return The asset ids removed, or the ones that would be, invisibly.
#' @examples
#' \dontrun{
#' gee_clean_assets(dry_run = TRUE)   # list what would be removed
#' }
#' @export
gee_clean_assets <- function(older_than_hours = 48, dry_run = FALSE, project = NULL) {
  ensure_gee_authenticated(project)
  ee <- reticulate::import("ee")
  project <- .resolve_project(project)
  assets <- tryCatch(ee$data$listAssets(list(parent = sprintf("projects/%s/assets", project)))[["assets"]],
                     error = function(e) NULL)
  ids  <- vapply(assets, function(a) a[["id"]], character(1))
  mine <- grep("alphasdm_[a-z]+_[0-9]{14}_", ids, value = TRUE)
  if (length(mine) == 0) {
    sdm_info("No leftover AlphaSDM assets.")
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
  drop  <- !is.na(age_h) & age_h > older_than_hours & !basename(mine) %in% active
  targets <- mine[drop]

  if (length(targets) == 0) {
    sdm_info(sprintf("%d AlphaSDM asset%s present, none old enough to remove.",
                     length(mine), if (length(mine) == 1) "" else "s"))
    return(invisible(character(0)))
  }
  if (!dry_run) ee_delete_asset_quietly(targets)
  sdm_info(sprintf("%s %d leftover AlphaSDM asset%s (older than %g h).",
                   if (dry_run) "Would remove" else "Removed", length(targets),
                   if (length(targets) == 1) "" else "s", older_than_hours))
  invisible(targets)
}
