#' Extract Alpha Earth Embeddings to Google Drive
#'
#' This script provides a powerful way to extract 64-band embeddings for a large number of coordinates.
#' It handles large datasets by splitting the export into multiple GEE tasks (one per year),
#' preventing "Out of Memory" (Error 8) failures on the GEE servers.
#'
#' @param df Data frame with columns: longitude, latitude, year
#' @param output_name Filename prefix for the output CSVs (year will be appended)
#' @param scale Resolution in meters (default 10)
#' @param drive_folder Google Drive folder name (default "AlphaSDM_embeddings")
#' @param wait Logical. If TRUE, waits for each task to queue before starting the next (recommended for big jobs).
#' @return A list of GEE Task objects
#'
#' @export
extract_embeddings_to_drive <- function(df, output_name, scale = 10, drive_folder = "AlphaSDM_embeddings", wait = TRUE) {
  # Direct source of internal logic
  if (file.exists("R/utils.R")) source("R/utils.R")
  if (file.exists("R/gee_logic.R")) source("R/gee_logic.R")
  library(rgee)
  ee <- reticulate::import("ee")

  # 1. Clean and deduplicate
  orig_n <- nrow(df)
  df <- df[!is.na(df$longitude) & !is.na(df$latitude) & !is.na(df$year), ]
  df <- unique(df[, c("longitude", "latitude", "year")])
  new_n <- nrow(df)
  
  message(sprintf("AlphaSDM: Extracting embeddings for %d unique coordinates (from %d original rows)...", new_n, orig_n))

  # 2. Get unique years
  years <- sort(unique(as.numeric(df$year)))
  message(sprintf("  -> Detected %d unique year(s): %s", length(years), paste(years, collapse = ", ")))
  
  tasks <- list()
  
  # 3. Process Year-by-Year to prevent OOM
  for (yr in years) {
    timestamp_message(sprintf("\n--- Processing Year: %d ---", yr))
    yr_df <- df[df$year == yr, ]
    
    # Upload points for this year only
    # Split into 50k chunks for safety
    CHUNK_SIZE <- 50000
    chunk_starts <- seq(1, nrow(yr_df), by = CHUNK_SIZE)
    yr_fcs <- list()
    for (i in seq_along(chunk_starts)) {
      start <- chunk_starts[i]
      end <- min(start + CHUNK_SIZE - 1, nrow(yr_df))
      message(sprintf("  -> Uploading batch %d/%d (%d points)...", i, length(chunk_starts), (end-start+1)))
      # We use our point-compression helper but we need to ensure unique IDs across chunks if needed
      # Actually for Drive export, geometries/IDs are optional if we keep lat/lon
      yr_fcs[[i]] <- upload_points_to_gee(yr_df[start:end, ], id_offset = start-1)
    }
    yr_full_fc <- ee$FeatureCollection(unname(yr_fcs))$flatten()
    
    # Sample
    yr_img <- get_embedding_image(yr, scale)
    sampled <- yr_img$sampleRegions(
      collection = yr_full_fc,
      properties = as.list(c("longitude", "latitude", "year")),
      scale = scale,
      tileScale = 16L, 
      geometries = FALSE
    )$filter(ee$Filter$notNull(as.list("A00")))

    # Export with Retry Logic (Ensures robustness against GEE server hiccups)
    task_name <- sprintf("%s_%d", output_name, yr)
    message(sprintf("  -> Queuing GEE Export: '%s'...", task_name))
    
    max_queue_retries <- 3
    for (attempt in seq_len(max_queue_retries)) {
      status <- try({
        task <- ee_table_to_drive(
          collection = sampled,
          description = task_name,
          folder = drive_folder,
          fileFormat = "CSV",
          selectors = as.list(c("longitude", "latitude", "year", sprintf("A%02d", 0:63)))
        )
        task$start()
      }, silent = TRUE)
      
      if (!inherits(status, "try-error")) {
        tasks[[as.character(yr)]] <- task
        break
      } else if (attempt == max_queue_retries) {
        stop("AlphaSDM: Failed to queue GEE task after ", max_queue_retries, " retry attempts. Error: ", as.character(status))
      }
      
      message("     ! GEE server hiccup. Retrying queue (Attempt ", attempt, "/", max_queue_retries, ")...")
      Sys.sleep(10) # 10s wait for server to breathe
    }
    
    if (wait) {
      message("  -> Task queued. Waiting for system handoff...")
      Sys.sleep(5) 
    }
  }

  message("\nSuccess! All Year-specific Export tasks are QUEUED on GEE.")
  message(sprintf("Total tasks: %d", length(tasks)))
  message("You can monitor progress with: rgee::ee_monitoring(task)")
  
  return(tasks)
}
