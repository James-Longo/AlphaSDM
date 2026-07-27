#' Format Data for AlphaSDM
#'
#' This function standardizes input data frames for use in the AlphaSDM pipeline.
#' It renames coordinate/year/presence columns to standard lowercase names,
#' preserves Alpha Earth embeddings if present, and removes all other columns.
#'
#' @param data A data frame containing your survey data.
#' @param coords A character vector of length 2 specifying the longitude and latitude columns IN ORDER: c(longitude_col, latitude_col). Note: Longitude first, then Latitude!
#' @param year A character string specifying the year or date column.
#' @param presence Optional. A character string specifying the presence/absence column (values should be 0 or 1).
#' @param species Optional. A character string specifying the species name column.
#' @return A standardized data frame ready for `AlphaSDM()` with lowercase column names.
#' @export
format_data <- function(data, coords, year, presence = NULL, species = NULL, label = NULL) {
    if (!isTRUE(.alphasdm_env$standardization_active)) {
        sdm_section("Data Standardization")
        .alphasdm_env$standardization_active <- TRUE
    }
    
    if (!is.data.frame(data)) {
        stop("Input 'data' must be a data frame.")
    }

    # 1. Coordinate Validation
    if (length(coords) != 2) {
        stop("'coords' must be a character vector of length 2: c(longitude_col, latitude_col)")
    }

    missing_cols <- setdiff(coords, names(data))
    if (length(missing_cols) > 0) {
        stop(paste("Coordinate columns not found in data:", paste(missing_cols, collapse = ", ")))
    }
    
    show_info <- !isTRUE(.alphasdm_env$standardization_info_printed)
    if (show_info) {
        sdm_info(sprintf("Validating coordinates: %s, %s", coords[1], coords[2]), indent = 1L)
    }

    # 2. Year Validation
    if (!year %in% names(data)) {
        stop(paste("Year column not found:", year))
    }

    # 3. Presence Validation
    if (!is.null(presence) && !presence %in% names(data)) {
        stop(paste("Presence column not found:", presence))
    }

    # 4. Species Validation
    if (!is.null(species) && !species %in% names(data)) {
        stop(paste("Species column not found:", species))
    }

    # 5. Build the output data frame with only required columns
    result <- data.frame(
        longitude = data[[coords[1]]],
        latitude = data[[coords[2]]],
        stringsAsFactors = FALSE
    )

    # Sanity check: Warn if coordinates appear swapped
    sample_lat <- result$latitude[!is.na(result$latitude)][1]
    sample_lon <- result$longitude[!is.na(result$longitude)][1]

    if (!is.null(sample_lat) && !is.null(sample_lon)) {
        lat_in_range <- sample_lat >= -90 && sample_lat <= 90
        lon_in_range <- sample_lon >= -180 && sample_lon <= 180

        if (!lat_in_range && lon_in_range) {
            warning(sprintf(
                "Coordinates may be SWAPPED! Latitude=%s is outside valid range [-90, 90].\n  Did you pass coords in the correct order? It should be: coords = c(longitude_col, latitude_col)\n  Your call: coords = c('%s', '%s')",
                sample_lat, coords[1], coords[2]
            ))
        }
    }

    # 6. Process Year (handle dates)
    year_data <- data[[year]]
    if (inherits(year_data, c("Date", "POSIXt"))) {
        result$year <- as.numeric(format(year_data, "%Y"))
    } else {
        # Try numeric conversion
        val <- suppressWarnings(as.numeric(year_data))

        if (all(is.na(val) & !is.na(year_data))) {
            # Likely date strings - try common formats
            formats <- c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d")
            parsed <- as.Date(rep(NA, nrow(data)))

            for (fmt in formats) {
                try_date <- as.Date(year_data, format = fmt)
                if (any(!is.na(try_date))) {
                    parsed <- try_date
                    break
                }
            }

            # Fallback to standard as.Date
            if (all(is.na(parsed))) {
                try_date <- try(as.Date(year_data), silent = TRUE)
                if (!inherits(try_date, "try-error")) {
                    parsed <- try_date
                }
            }

            result$year <- as.numeric(format(parsed, "%Y"))
        } else {
            result$year <- val
        }
    }
    if (show_info) {
        sdm_info(sprintf("Standardizing time/dates using column: %s", year), indent = 1L)
    }

    # 8. Add presence column (lowercase "present")
    if (!is.null(presence)) {
        result$present <- as.numeric(data[[presence]])
    } else {
        sdm_warn("No 'presence' column supplied; treating all records as presence-only")
        sdm_info("Background pseudo-absences will be generated automatically", indent = 1L)
        result$present <- 1
    }

    # 9. Add species column (lowercase "species")
    if (!is.null(species)) {
        result$species <- as.character(data[[species]])
    }

    # 11. (REMOVED) Local preservation of Alpha Earth Embeddings is disabled per security policy.
    # Prediction and training always happen server-side on GEE.

    # 12. Filter to years with Alpha Earth data (2017+)
    if (show_info) {
        sdm_info("Filtering for Alpha Earth coverage (2017\u20132025)", indent = 1L)
        .alphasdm_env$standardization_info_printed <- TRUE
    }
    rows_before <- nrow(result)
    result <- result[result$year >= 2017 & result$year <= 2025, ]
    rows_after <- nrow(result)

    if (rows_before != rows_after) {
        sdm_warn(sprintf("%d row%s removed \u2014 outside Alpha Earth temporal coverage (2017\u20132025)",
                         rows_before - rows_after,
                         if ((rows_before - rows_after) == 1) "" else "s"))
    }

    # 13. Remove rows with NA values
    rows_before <- nrow(result)
    result <- na.omit(result)
    rows_after <- nrow(result)

    if (rows_before != rows_after) {
        sdm_warn(sprintf("%d row%s removed \u2014 contains missing (NA) values",
                         rows_before - rows_after,
                         if ((rows_before - rows_after) == 1) "" else "s"))
    }

    # 14. Remove duplicate rows (Prioritize presence: if coords/year collide, take max(present))
    #
    # The presence-first sort is also the pipeline's row-order contract: GEE's libsvm
    # assigns its positive class from the FIRST label it sees in the training data, and
    # predict_gee_map's SVM flip assumes presence (1) is positive. Emitting presences
    # first here keeps that invariant at the data boundary, so no downstream stage needs
    # to re-sort. (order() is stable, so original order is preserved within each class.)
    rows_before <- nrow(result)

    key_cols <- setdiff(names(result), "present")
    result <- result[order(-result$present), ]
    result <- result[!duplicated(result[, key_cols, drop = FALSE]), ]
    rownames(result) <- NULL
    rows_after <- nrow(result)

    if (rows_before != rows_after) {
        sdm_warn(sprintf("%d row%s removed \u2014 duplicate or conflicting records (presence prioritized)",
                         rows_before - rows_after,
                         if ((rows_before - rows_after) == 1) "" else "s"))
    }

    # List final columns
    cols_desc <- paste(names(result), collapse = ", ")
    prefix <- if (!is.null(label)) paste0(label, " data") else "Data"
    sdm_done(sprintf("%s ready: %d rows × %d columns (%s)", prefix, nrow(result), ncol(result), cols_desc))

    return(result)
}
