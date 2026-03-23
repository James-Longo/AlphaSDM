#' Extract Alpha Earth Embeddings to Google Drive
#'
#' This script provides a powerful way to extract 64-band embeddings for a large number of coordinates.
#' It deduplicates coordinates, uploads them to GEE in compressed batches, and starts a Drive export task.
#'
#' @param df Data frame with columns: longitude, latitude, year
#' @param output_name Filename for the output CSV (without extension)
#' @param scale Resolution in meters (default 10)
#' @param drive_folder Google Drive folder name (default "embedding output")
#' @return GEE Task object
#'
#' @examples
#' \dontrun{
#'   unique_points <- PO_train %>% distinct(lon, lat, year) %>% rename(longitude=lon, latitude=lat)
#'   task <- extract_embeddings_to_drive(unique_points, "glc25_po_embeddings")
#'   # Monitor with: rgee::ee_monitoring(task)
#' }
extract_embeddings_to_drive <- function(df, output_name, scale = 10, drive_folder = "embedding output") {
  # Direct source of internal logic for script robustness
  if (file.exists("R/utils.R")) source("R/utils.R")
  if (file.exists("R/gee_logic.R")) source("R/gee_logic.R")
  library(rgee)
  ee <- reticulate::import("ee")

  # 1. Clean and deduplicate (essential for performance)
  orig_n <- nrow(df)
  df <- df[!is.na(df$longitude) & !is.na(df$latitude) & !is.na(df$year), ]
  df <- unique(df[, c("longitude", "latitude", "year")])
  new_n <- nrow(df)
  
  message(sprintf("Extracting embeddings for %d unique coordinates (from %d original rows)...", new_n, orig_n))

  # 2. Upload points to GEE (Uses autoSDM's high-efficiency MultiPoint compression)
  # We split into 100k chunks to avoid getting hit by GEE's 10MB payload limit
  MAX_UPLOAD_CHUNK <- 100000 
  chunk_starts <- seq(1, nrow(df), by = MAX_UPLOAD_CHUNK)
  
  all_fcs <- list()
  for (i in seq_along(chunk_starts)) {
    start <- chunk_starts[i]
    end <- min(start + MAX_UPLOAD_CHUNK - 1, nrow(df))
    chunk_df <- df[start:end, , drop = FALSE]
    
    message(sprintf("  -> Uploading batch %d/%d (%d points)...", i, length(chunk_starts), nrow(chunk_df)))
    all_fcs[[i]] <- upload_points_to_gee(chunk_df, id_offset = start - 1)
  }
  
  full_fc <- ee$FeatureCollection(unname(all_fcs))$flatten()

  # 3. Sample Embeddings
  # We use the internal autoSDM function which handles annual image mosaicking
  message("  -> Sampling embeddings from Alpha Earth (GEE-side)...")
  sampled_fc <- get_embeddings_at_fc(
    fc = full_fc, 
    scale = scale, 
    properties = c("year", "longitude", "latitude"),
    geometries = FALSE
  )

  # 4. Start Export Task
  message(sprintf("  -> Starting GEE Export task: '%s' to folder: '%s'...", output_name, drive_folder))
  task <- ee_table_to_drive(
    collection = sampled_fc,
    description = output_name,
    folder = drive_folder,
    fileFormat = "CSV",
    selectors = as.list(c("longitude", "latitude", "year", sprintf("A%02d", 0:63)))
  )
  
  task$start()
  
  message("\nSuccess! Export task is QUEUED on GEE.")
  message("You can monitor progress with: rgee::ee_monitoring(task)")
  message("Or check your GEE Task tab at: https://code.earthengine.google.com/")
  
  return(task)
}

# --- Quick Test Setup ---
# if (interactive()) {
#   ee_activate()
#   # test_df <- data.frame(longitude=c(2.3522, 13.4050), latitude=c(48.8566, 52.5200), year=c(2021, 2021))
#   # extract_embeddings_to_drive(test_df, "autoSDM_test_extract")
# }
