# Access to the concept label table and visual verification of concepts.
# The table is produced once per SAE version by dev/03_autolabel.R and ships
# as package data; the raw correlate columns in it are the ground truth and
# the label column is a convenience.

#' Read the shipped concept label table
#'
#' One row per SAE latent: `concept_id`, `label`, the top catalog correlates
#' behind that label, a `coherence` score (how strongly the best correlate
#' explains the latent), and `example_locations`, the top-activating points as
#' `"lon,lat;lon,lat;..."`.
#'
#' @param path Path to a label CSV written by dev/03_autolabel.R; NULL for the
#'   table installed with the package.
#' @return A data frame with one row per concept.
#' @export
concept_labels <- function(path = NULL) {
  if (is.null(path)) {
    path <- system.file("extdata", "concept_labels_v1.csv", package = "AlphaSDM")
    if (!nzchar(path)) {
      stop("No concept label table is installed with this copy of AlphaSDM. ",
           "The concept layer is experimental: build the table with ",
           "dev/03_autolabel.R, or pass `path = `.", call. = FALSE)
    }
  }
  labels <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!"concept_id" %in% names(labels))
    stop("The label table has no concept_id column: ", path, call. = FALSE)
  labels
}

#' Parse an example_locations string into a lon/lat data frame
#' @noRd
parse_example_locations <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(NULL)
  pairs <- strsplit(strsplit(s, ";", fixed = TRUE)[[1]], ",", fixed = TRUE)
  df <- do.call(rbind, lapply(pairs, function(p) {
    v <- suppressWarnings(as.numeric(p))
    if (length(v) != 2L || anyNA(v)) return(NULL)
    data.frame(longitude = v[1], latitude = v[2])
  }))
  if (is.null(df) || nrow(df) == 0L) NULL else df
}

#' Show Sentinel-2 chips at a concept's top-activating locations
#'
#' Pulls RGB thumbnails at the top-activating points recorded in the label
#' table, for eyeball verification that a concept's label matches what is on
#' the ground.
#'
#' @param concept_id Concept id, e.g. `"C017"`.
#' @param n How many chips to show, up to the number on record.
#' @param year Year of Sentinel-2 imagery to composite.
#' @param buffer_m Chip half-width in metres around each point.
#' @param thumb_px Thumbnail edge length in pixels.
#' @param max_cloud Scenes above this CLOUDY_PIXEL_PERCENTAGE are excluded
#'   from the composite.
#' @param viz_max Reflectance mapped to full brightness (standard Sentinel-2
#'   true-colour stretch).
#' @param labels Optional label table; default the one installed with the package.
#' @return The chip file paths, invisibly. Draws the chips as a side effect.
#' @export
plot_concept_chips <- function(concept_id, n = 9L, year = 2023,
                               buffer_m = 500, thumb_px = 256L,
                               max_cloud = 20, viz_max = 3000,
                               labels = NULL) {
  if (!requireNamespace("png", quietly = TRUE))
    stop("plot_concept_chips() needs the 'png' package: install.packages('png')",
         call. = FALSE)
  ensure_gee_authenticated()
  ee <- reticulate::import("ee")

  if (is.null(labels)) labels <- concept_labels()
  row <- labels[labels$concept_id == concept_id, , drop = FALSE]
  if (nrow(row) == 0L)
    stop("Concept '", concept_id, "' is not in the label table.", call. = FALSE)
  locs <- parse_example_locations(row$example_locations[1L])
  if (is.null(locs))
    stop("Concept '", concept_id, "' has no example locations on record.", call. = FALSE)
  locs <- locs[seq_len(min(as.integer(n), nrow(locs))), , drop = FALSE]

  s2 <- ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
    filterDate(sprintf("%d-01-01", year), sprintf("%d-01-01", year + 1))$
    filter(ee$Filter$lt("CLOUDY_PIXEL_PERCENTAGE", max_cloud))

  paths <- character(0)
  for (i in seq_len(nrow(locs))) {
    region <- ee$Geometry$Point(c(locs$longitude[i], locs$latitude[i]))$
      buffer(buffer_m)$bounds()
    img <- s2$filterBounds(region)$median()$select(c("B4", "B3", "B2"))
    url <- retry_curl_download(img$getThumbURL(list(
      region = region, dimensions = as.integer(thumb_px),
      min = 0, max = viz_max, format = "png")))
    dest <- tempfile(fileext = ".png")
    ok <- tryCatch({ utils::download.file(url, dest, mode = "wb", quiet = TRUE); TRUE },
                   error = function(e) FALSE)
    if (ok) paths <- c(paths, dest)
  }
  if (length(paths) == 0L)
    stop("No chips could be downloaded for '", concept_id, "'.", call. = FALSE)

  side <- ceiling(sqrt(length(paths)))
  op <- graphics::par(mfrow = c(side, side), mar = c(0.2, 0.2, 1.2, 0.2))
  on.exit(graphics::par(op), add = TRUE)
  for (i in seq_along(paths)) {
    im <- png::readPNG(paths[i])
    graphics::plot.new()
    graphics::rasterImage(im, 0, 0, 1, 1)
    graphics::title(sprintf("%s #%d", concept_id, i), cex.main = 0.8)
  }
  invisible(paths)
}
