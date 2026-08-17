# Google Drive access for image exports. toDrive is the only export destination
# with shardSize (compute memory) and fileDimensions (output tiling). No separate
# sign-in: the Earth Engine credential already carries the drive scope.

.drive <- new.env(parent = emptyenv())

#' Authorized HTTP session for the Drive API
#'
#' @return A `google.auth.transport.requests.AuthorizedSession`, cached per session.
#' @noRd
drive_session <- function() {
  if (!is.null(.drive$session)) return(.drive$session)
  ee    <- reticulate::import("ee")
  creds <- reticulate::import("google.oauth2.credentials")
  treq  <- reticulate::import("google.auth.transport.requests")
  json  <- reticulate::import("json")
  bi    <- reticulate::import_builtins()

  path <- ee$oauth$get_credentials_path()
  if (!file.exists(path))
    stop("No Earth Engine credentials found. Run setup_gee() first.", call. = FALSE)
  fh  <- bi$open(path); cfg <- json$load(fh); fh$close()
  if (is.null(cfg[["refresh_token"]]))
    stop("The stored Earth Engine credential has no refresh token. Re-run setup_gee().",
         call. = FALSE)

  cr <- creds$Credentials(
    NULL, refresh_token = cfg[["refresh_token"]],
    token_uri = "https://oauth2.googleapis.com/token",
    client_id = ee$oauth$CLIENT_ID, client_secret = ee$oauth$CLIENT_SECRET,
    scopes = cfg[["scopes"]])
  cr$refresh(treq$Request())
  .drive$session <- treq$AuthorizedSession(cr)
  .drive$session
}

#' List the files an export has written so far
#'
#' @param prefix File name prefix, which is the export's `fileNamePrefix`.
#' @return A data frame with `id`, `name` and `bytes`, empty when nothing is written.
#' @noRd
drive_list_files <- function(prefix) {
  s <- drive_session()
  # Drive's "name contains" silently returns nothing for some long compound
  # strings, so query on a short leading token and filter exactly here.
  q_token <- paste(head(strsplit(prefix, "_", fixed = TRUE)[[1]], 2L), collapse = "_")
  out <- list(); token <- NULL
  repeat {
    params <- list(q = sprintf("name contains '%s' and trashed = false", q_token),
                   fields = "nextPageToken,files(id,name,size)",
                   pageSize = 1000L)
    if (!is.null(token)) params$pageToken <- token
    resp <- s$get("https://www.googleapis.com/drive/v3/files", params = params)
    d <- resp$json()
    fl <- d[["files"]]
    if (!is.null(fl)) out <- c(out, fl)
    token <- d[["nextPageToken"]]
    if (is.null(token)) break
  }
  if (!length(out))
    return(data.frame(id = character(0), name = character(0), bytes = numeric(0),
                      stringsAsFactors = FALSE))
  df <- data.frame(
    id    = vapply(out, function(f) as.character(f[["id"]]), character(1)),
    name  = vapply(out, function(f) as.character(f[["name"]]), character(1)),
    bytes = vapply(out, function(f) {
      v <- suppressWarnings(as.numeric(f[["size"]])); if (length(v) != 1L) NA_real_ else v
    }, numeric(1)),
    stringsAsFactors = FALSE)
  df[startsWith(df$name, prefix), , drop = FALSE]
}

#' Download one Drive file
#'
#' The body is written from Python so the bytes never cross into R, which matters
#' for tiles that run to tens of megabytes.
#'
#' @noRd
drive_download_file <- function(id, dest) {
  s  <- drive_session()
  bi <- reticulate::import_builtins()
  resp <- s$get(sprintf("https://www.googleapis.com/drive/v3/files/%s", id),
                params = list(alt = "media"))
  if (as.integer(resp$status_code) != 200L)
    stop(sprintf("Drive refused the download of %s (HTTP %s)", id, resp$status_code),
         call. = FALSE)
  fh <- bi$open(dest, "wb"); fh$write(resp$content); fh$close()
  dest
}

#' Decide whether a Drive file name was written by a given export
#'
#' Matched files are deleted after download, and Drive's "name contains" query is
#' a loose substring match, so accept only names Earth Engine writes: the prefix,
#' an optional _rNNN band suffix, then the -yMin-xMin tile suffix or nothing for a
#' single-file export.
#'
#' @param prefix The export's fileNamePrefix.
#' @param names Character vector of Drive file names.
#' @return Logical vector, TRUE for names the export wrote.
#' @noRd
is_export_tile_name <- function(prefix, names) {
  grepl(sprintf("^%s(_r\\d+(_c\\d+)?)?(-\\d+-\\d+)?\\.tif$", prefix), names)
}

#' Delete one Drive file
#' @noRd
drive_delete_file <- function(id) {
  s <- drive_session()
  invisible(tryCatch(s$delete(sprintf("https://www.googleapis.com/drive/v3/files/%s", id)),
                     error = function(e) NULL))
}

#' Compute an image through Drive and bring the tiles back
#'
#' Earth Engine computes the region in the batch system, splits the result into
#' aligned tiles of `file_dim` pixels, and writes them to Drive. The tiles are then
#' downloaded and, unless `keep_on_drive` is TRUE, removed from Drive so a large map
#' does not sit in the user's storage.
#'
#' The image is left masked rather than filled with a sentinel, so `skipEmptyTiles`
#' can drop tiles that are entirely masked. On a coastal or island region that is
#' most of the bounding box.
#'
#' @param image  ee.Image with a single prediction band.
#' @param region ee.Geometry whose bounds define the extent.
#' @param scale  Output pixel size in metres.
#' @param out_dir Directory to write the tiles into. Created if absent.
#' @param folder Drive folder name.
#' @param shard_size Pixel block Earth Engine computes in. Its default is 256.
#' @param file_dim Pixel size of each output tile. Must be a multiple of `shard_size`.
#' @param poll_seconds Seconds between progress reports.
#' @param keep_on_drive TRUE leaves the tiles in Drive after downloading them.
#' @return `out_dir`, containing one GeoTIFF per tile.
#' @noRd
ee_export_image_drive <- function(image, region, scale, out_dir,
                                  folder = "AlphaSDM", shard_size = 256L,
                                  file_dim = 4096L, poll_seconds = 15,
                                  keep_on_drive = FALSE, valid_mask = NULL) {
  ee <- reticulate::import("ee")
  if (file_dim %% shard_size != 0L)
    stop(sprintf("file_dim (%d) must be a multiple of shard_size (%d).", file_dim,
                 shard_size), call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- sprintf("alphasdm_%s_%s", format(Sys.time(), "%Y%m%d%H%M%S"),
                    paste(sample(c(letters, 0:9), 6, replace = TRUE), collapse = ""))

  # One export task per square grid cell of at most two output files per side.
  # Small independent tasks keep a restart or a terminal failure cheap and give
  # each cell its own retry budget; a shared crsTransform keeps every cell on one
  # global pixel grid.
  ring <- region$bounds()$coordinates()$get(0L)$getInfo()
  xs <- vapply(ring, function(p) p[[1]], numeric(1))
  ys <- vapply(ring, function(p) p[[2]], numeric(1))
  west <- min(xs); east <- max(xs); south <- min(ys); north <- max(ys)
  dpp <- scale / 111320                       # degrees per pixel, as build_grid uses
  cell_px  <- 2L * file_dim
  cell_deg <- cell_px * dpp
  n_rows <- max(1L, as.integer(ceiling((north - south) / cell_deg)))
  n_cols <- max(1L, as.integer(ceiling((east - west) / cell_deg)))
  # ALPHASDM_BAND_GEOM: "transform" (default; shared grid) or "scale".
  # ALPHASDM_BAND_MODE: "rows" (default grid) or "single" for one whole-region task.
  use_transform <- !identical(Sys.getenv("ALPHASDM_BAND_GEOM", "transform"), "scale")
  if (identical(Sys.getenv("ALPHASDM_BAND_MODE", "rows"), "single")) { n_rows <- 1L; n_cols <- 1L }
  single <- n_rows == 1L && n_cols == 1L
  transform <- list(dpp, 0, west, 0, -dpp, north)

  # Submit only cells that intersect the caller's region geometry. The region is
  # the user's own statement of where the map should exist, so cells outside it
  # cost nothing; all-masked cells inside it are reclassified as empty on
  # completion.
  has_land <- matrix(TRUE, n_rows, n_cols)
  if (n_rows * n_cols > 1L) {
    feats <- list()
    for (b in seq_len(n_rows)) for (cc in seq_len(n_cols)) {
      top  <- north - (b - 1L) * cell_deg; bot <- max(south, top - cell_deg)
      left <- west + (cc - 1L) * cell_deg; rgt <- min(east, left + cell_deg)
      feats[[length(feats) + 1L]] <- ee$Feature(
        ee$Geometry$Rectangle(c(left, bot, rgt, top)), list(b = b, cc = cc))
    }
    hits <- tryCatch(
      ee$FeatureCollection(feats)$filterBounds(region)$
        reduceColumns(ee$Reducer$toList(2L), list("b", "cc"))$get("list")$getInfo(),
      error = function(e) NULL)
    if (!is.null(hits)) {
      has_land <- matrix(FALSE, n_rows, n_cols)
      for (hh in hits) has_land[hh[[1]], hh[[2]]] <- TRUE
      sdm_info(sprintf("%d of %d cells intersect the AOI; the rest are skipped.",
                       sum(has_land), n_rows * n_cols), indent = 2L)
    }
  }

  tasks <- list()
  for (b in seq_len(n_rows)) for (cc in seq_len(n_cols)) {
    if (!has_land[b, cc]) next
    top  <- if (single) north else north - (b - 1L) * cell_deg
    bot  <- if (single) south else max(south, top - cell_deg)
    left <- if (single) west  else west + (cc - 1L) * cell_deg
    rgt  <- if (single) east  else min(east, left + cell_deg)
    cell_prefix <- if (single) prefix else sprintf("%s_r%03d_c%03d", prefix, b, cc)
    task <- ee$batch$Export$image$toDrive(
      image = image, description = cell_prefix, folder = folder,
      fileNamePrefix = cell_prefix,
      region = ee$Geometry$Rectangle(c(left, bot, rgt, top)),
      crs = "EPSG:4326",
      crsTransform = if (use_transform) transform else NULL,
      scale = if (use_transform) NULL else scale,
      shardSize = shard_size, fileDimensions = file_dim, skipEmptyTiles = TRUE,
      maxPixels = 1e13, fileFormat = "GeoTIFF")
    task$start()
    tasks[[cell_prefix]] <- task
  }
  sdm_info(sprintf(
    "Computing the map on Earth Engine as %d task%s, writing tiles to Drive/%s",
    length(tasks), if (length(tasks) == 1L) "" else "s", folder), indent = 1L)
  ee_await_export_all(tasks, poll_seconds)

  files <- drive_list_files(prefix)
  files <- files[is_export_tile_name(prefix, files$name), , drop = FALSE]
  # A band of open water legitimately writes nothing; only a map with no tiles at
  # all is an error.
  if (!nrow(files))
    stop("The export finished but wrote no tiles. If the whole region is masked, ",
         "skipEmptyTiles will have dropped every tile.", call. = FALSE)

  total <- sum(files$bytes, na.rm = TRUE)
  sdm_info(sprintf("Downloading %d tile%s (%s)", nrow(files),
                   if (nrow(files) == 1L) "" else "s",
                   if (total >= 1024^3) sprintf("%.1f GB", total / 1024^3)
                   else sprintf("%.0f MB", total / 1024^2)), indent = 1L)
  for (i in seq_len(nrow(files))) {
    dest <- file.path(out_dir, files$name[i])
    drive_download_file(files$id[i], dest)
    got <- file.info(dest)$size
    # Compare against the size Drive reported. A short file means a truncated
    # download, and deleting the Drive copy would make it unrecoverable.
    if (!is.na(files$bytes[i]) && !is.na(got) && got != files$bytes[i])
      stop(sprintf("Tile %s downloaded as %.0f bytes but Drive reports %.0f.",
                   files$name[i], got, files$bytes[i]), call. = FALSE)
    if (!keep_on_drive) drive_delete_file(files$id[i])
    if (i %% 25 == 0 || i == nrow(files))
      sdm_info(sprintf("%d of %d tile%s downloaded", i, nrow(files),
                       if (nrow(files) == 1L) "" else "s"), indent = 2L)
  }
  sdm_info(sprintf("Tiles written to %s", out_dir), indent = 1L)
  out_dir
}

#' Wait for a set of export tasks, reporting joint progress
#'
#' One poll loop over every task; all scheduling is Earth Engine's. Waits for all
#' tasks to reach a terminal state before raising any failure, so completed bands
#' keep their output.
#'
#' @param tasks Named list of started tasks, name = description.
#' @param poll_seconds Seconds between polls.
#' @return Invisibly TRUE. Errors if any task failed or was cancelled.
#' @noRd
ee_await_export_all <- function(tasks, poll_seconds = 15) {
  ee <- reticulate::import("ee")
  start <- Sys.time(); last_beat <- -Inf
  # One listOperations call per poll covers every task; per-task status calls
  # would multiply requests by the task count.
  states_by_desc <- function() {
    out <- list()
    ops <- tryCatch(ee$data$listOperations(), error = function(e) list())
    for (o in ops) {
      m <- o$metadata
      d <- tryCatch(as.character(m$description), error = function(e) NULL)
      st <- tryCatch(as.character(m$state), error = function(e) NULL)
      if (length(d) == 1L && length(st) == 1L) out[[d]] <- st
    }
    out
  }
  ee_to_task_state <- c(PENDING = "READY", RUNNING = "RUNNING",
                        SUCCEEDED = "COMPLETED", COMPLETED = "COMPLETED",
                        FAILED = "FAILED", CANCELLED = "CANCELLED",
                        CANCELLING = "CANCELLED")
  repeat {
    sm <- states_by_desc()
    sts <- vapply(names(tasks), function(nm) {
      st <- sm[[nm]]
      if (is.null(st)) tasks[[nm]]$status()[["state"]]
      else { mapped <- ee_to_task_state[[st]]; if (is.null(mapped)) st else mapped }
    }, character(1))
    done <- sum(sts == "COMPLETED")
    bad  <- sum(sts %in% c("FAILED", "CANCELLED"))
    if (done + bad == length(sts)) break
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    if (elapsed - last_beat >= 60) {
      run <- which(sts == "RUNNING")
      detail <- if (length(run))
        ee_task_progress(tasks[[run[1]]]$status()[["name"]]) else ""
      # The bracket describes one running task; label it so its percent is not
      # read as the whole job's.
      if (nzchar(detail)) detail <- sub(" \\[", " [current task: ", detail)
      sdm_task_monitor_hint()
      sdm_info(sprintf("%d of %d tasks done, %d running, %d queued%s%s (%s elapsed)",
                       done, length(sts), length(run),
                       sum(sts == "READY"),
                       if (bad > 0L) sprintf(", %d FAILED", bad) else "", detail,
                       if (elapsed < 600) sprintf("%.0fs", elapsed)
                       else sprintf("%.0f min", elapsed / 60)), indent = 2L)
      last_beat <- elapsed
    }
    Sys.sleep(poll_seconds)
  }
  if (bad > 0L) {
    who <- names(tasks)[sts %in% c("FAILED", "CANCELLED")]
    msgs <- vapply(who, function(nm) {
      s <- tasks[[nm]]$status()
      if (!is.null(s[["error_message"]])) s[["error_message"]] else s[["state"]]
    }, character(1))
    # A cell that is entirely masked (open ocean) is reported FAILED by Earth
    # Engine but is an empty result, not a failure.
    empty <- grepl("No valid \\(un-masked\\) pixels", msgs)
    if (any(empty))
      sdm_info(sprintf("%d cell%s contained no unmasked pixels (open water); skipped.",
                       sum(empty), if (sum(empty) == 1L) "" else "s"), indent = 2L)
    who <- who[!empty]; msgs <- msgs[!empty]
    if (length(who) > 0L)
      stop(sprintf("%d of %d export tasks did not complete. %s. Completed cells remain in Drive.",
                   length(who), length(sts),
                   paste(sprintf("%s: %s", who, msgs), collapse = "; ")), call. = FALSE)
  }
  invisible(TRUE)
}
