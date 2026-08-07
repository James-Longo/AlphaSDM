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
  out <- list(); token <- NULL
  repeat {
    params <- list(q = sprintf("name contains '%s' and trashed = false", prefix),
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
  data.frame(
    id    = vapply(out, function(f) as.character(f[["id"]]), character(1)),
    name  = vapply(out, function(f) as.character(f[["name"]]), character(1)),
    bytes = vapply(out, function(f) {
      v <- suppressWarnings(as.numeric(f[["size"]])); if (length(v) != 1L) NA_real_ else v
    }, numeric(1)),
    stringsAsFactors = FALSE)
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
#' a loose substring match, so accept only names Earth Engine writes: the prefix
#' plus the -yMin-xMin tile suffix, or the bare prefix for a single-file export.
#'
#' @param prefix The export's fileNamePrefix.
#' @param names Character vector of Drive file names.
#' @return Logical vector, TRUE for names the export wrote.
#' @noRd
is_export_tile_name <- function(prefix, names) {
  grepl(sprintf("^%s(-\\d+-\\d+)?\\.tif$", prefix), names)
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
                                  keep_on_drive = FALSE) {
  ee <- reticulate::import("ee")
  if (file_dim %% shard_size != 0L)
    stop(sprintf("file_dim (%d) must be a multiple of shard_size (%d).", file_dim,
                 shard_size), call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  prefix <- sprintf("alphasdm_%s_%s", format(Sys.time(), "%Y%m%d%H%M%S"),
                    paste(sample(c(letters, 0:9), 6, replace = TRUE), collapse = ""))
  task <- ee$batch$Export$image$toDrive(
    image = image, description = prefix, folder = folder, fileNamePrefix = prefix,
    region = region, scale = scale, crs = "EPSG:4326",
    shardSize = shard_size, fileDimensions = file_dim, skipEmptyTiles = TRUE,
    maxPixels = 1e13, fileFormat = "GeoTIFF")
  task$start()
  sdm_info(sprintf("Computing the map on Earth Engine and writing tiles to Drive/%s",
                   folder), indent = 1L)
  ee_await_export(list(task = task, asset_id = prefix, drive_prefix = prefix),
                  poll_seconds)

  files <- drive_list_files(prefix)
  files <- files[is_export_tile_name(prefix, files$name), , drop = FALSE]
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
