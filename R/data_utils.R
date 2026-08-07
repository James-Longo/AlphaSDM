#' Standardize occurrence records for AlphaSDM
#'
#' Standardizes an input data frame for the AlphaSDM pipeline: renames the
#' coordinate, year and presence columns to the package's lowercase names, drops
#' everything else, and filters to the years Alpha Earth covers.
#'
#' The embeddings are annual. Records dated outside the covered window are dropped
#' and the number removed is reported; set them to the first covered year to keep
#' them instead. The window is read from the Earth Engine collection, so it tracks
#' each annual release. This means `format_data()` needs a connection: run
#' [setup_gee()] first.
#'
#' @param data A data frame of survey records.
#' Coordinates must be longitude and latitude in decimal degrees on WGS84
#' (EPSG:4326), the system GPS units and most occurrence databases (GBIF, eBird)
#' report in. A data frame carries no CRS information, so nothing can be
#' reprojected on your behalf: projected coordinates such as UTM metres are
#' rejected, and coordinates in another geographic datum are not detectable and
#' will be treated as WGS84. Reproject with
#' `sf::st_transform(x, 4326)` before formatting if needed.
#'
#' @param coords Character vector of length 2 naming the longitude and latitude
#'   columns, longitude first: `c(longitude_col, latitude_col)`. Values must be
#'   decimal degrees on WGS84 (EPSG:4326); see Details.
#' @param year A character string naming the year or date column. Dates are reduced
#'   to their year. Records outside the Alpha Earth window are dropped.
#' @param presence Optional. Name of the presence column, holding 1 for presence
#'   and 0 for absence. Omit it to treat every record as a presence.
#' @param species Optional. Name of the species column.
#' @param label Optional. Short name for this data set, used only to label the
#'   console summary line, for example "Training" or "Evaluation".
#' @return A data frame with `longitude`, `latitude`, `year` and `present` columns,
#'   ready for [evaluate_models()] or [generate_map()]. Rows are ordered presences
#'   first. Fewer rows come back than went in whenever records are dropped for
#'   coverage, missing values or duplication; each drop is reported as a message.
#' @export
format_data <- function(data, coords, year, presence = NULL, species = NULL, label = NULL) {
    if (!isTRUE(.alphasdm_env$standardization_active)) {
        sdm_section("Data Standardization")
        .alphasdm_env$standardization_active <- TRUE
    }
    
    if (!is.data.frame(data)) {
        stop("Input 'data' must be a data frame.")
    }

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

    if (!year %in% names(data)) {
        stop(paste("Year column not found:", year))
    }

    if (!is.null(presence) && !presence %in% names(data)) {
        stop(paste("Presence column not found:", presence))
    }

    if (!is.null(species) && !species %in% names(data)) {
        stop(paste("Species column not found:", species))
    }

    result <- data.frame(
        longitude = data[[coords[1]]],
        latitude = data[[coords[2]]],
        stringsAsFactors = FALSE
    )

    # Catch swapped coordinates: a latitude outside [-90, 90] beside a plausible
    # longitude almost always means the two columns were given the wrong way round.
    sample_lat <- result$latitude[!is.na(result$latitude)][1]
    sample_lon <- result$longitude[!is.na(result$longitude)][1]

    # Guard on NA, not NULL: `x[1]` on an empty vector gives NA. An all-NA
    # coordinate column would otherwise reach the comparison below and fail with
    # "missing value where TRUE/FALSE needed".
    # Out-of-range coordinates never proceed: Earth Engine would treat them as
    # degrees regardless. Latitude out of range beside a plausible longitude
    # usually means the columns were passed in the wrong order; out of range in
    # both usually means a projected system such as UTM metres.
    bad_lat <- !is.na(result$latitude)  & (result$latitude  < -90  | result$latitude  > 90)
    bad_lon <- !is.na(result$longitude) & (result$longitude < -180 | result$longitude > 180)
    if (any(bad_lat) || any(bad_lon)) {
        if (any(bad_lat) && !any(bad_lon)) {
            stop(sprintf(paste0(
                "Coordinates may be SWAPPED: latitude=%s is outside [-90, 90] but the ",
                "longitude column looks valid. coords must be longitude first: ",
                "coords = c(longitude_col, latitude_col). Your call: coords = c('%s', '%s')"),
                result$latitude[bad_lat][1], coords[1], coords[2]), call. = FALSE)
        }
        stop(sprintf(paste0(
            "%d record(s) have coordinates outside decimal-degree range ",
            "(longitude [-180, 180], latitude [-90, 90]). AlphaSDM needs longitude/latitude ",
            "in decimal degrees on WGS84 (EPSG:4326). If your coordinates are projected ",
            "(e.g. UTM metres), reproject first: sf::st_transform(x, 4326)."),
            sum(bad_lat | bad_lon)), call. = FALSE)
    }

    year_data <- data[[year]]
    if (inherits(year_data, c("Date", "POSIXt"))) {
        result$year <- as.numeric(format(year_data, "%Y"))
    } else {
        val <- suppressWarnings(as.numeric(year_data))

        if (all(is.na(val) & !is.na(year_data))) {
            # Nothing parsed as a number, so treat the column as date strings.
            formats <- c("%Y-%m-%d", "%m/%d/%Y", "%d/%m/%Y", "%Y/%m/%d")
            parsed <- as.Date(rep(NA, nrow(data)))

            for (fmt in formats) {
                try_date <- as.Date(year_data, format = fmt)
                if (any(!is.na(try_date))) {
                    parsed <- try_date
                    break
                }
            }

            # None of the explicit formats matched, so let as.Date guess.
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

    if (!is.null(presence)) {
        result$present <- as.numeric(data[[presence]])
    } else {
        sdm_warn("No 'presence' column supplied; treating all records as presence-only")
        sdm_info("Background pseudo-absences will be generated automatically", indent = 1L)
        result$present <- 1
    }

    if (!is.null(species)) {
        result$species <- as.character(data[[species]])
    }

    # Filter to the years Alpha Earth covers, read from the collection itself.
    yrs <- alphaearth_year_range()
    if (show_info) {
        sdm_info(sprintf("Filtering for Alpha Earth coverage (%d-%d)", yrs[1], yrs[2]), indent = 1L)
        .alphasdm_env$standardization_info_printed <- TRUE
    }
    rows_before <- nrow(result)
    result <- result[result$year >= yrs[1] & result$year <= yrs[2], ]
    rows_after <- nrow(result)

    if (rows_before != rows_after) {
        sdm_warn(sprintf("%d row%s removed: outside Alpha Earth temporal coverage (%d-%d). Set them to %d to keep them.",
                         rows_before - rows_after,
                         if ((rows_before - rows_after) == 1) "" else "s",
                         yrs[1], yrs[2], yrs[1]))
    }

    rows_before <- nrow(result)
    result <- na.omit(result)
    rows_after <- nrow(result)

    if (rows_before != rows_after) {
        sdm_warn(sprintf("%d row%s removed: contains missing (NA) values",
                         rows_before - rows_after,
                         if ((rows_before - rows_after) == 1) "" else "s"))
    }

    # Drop duplicates, keeping the presence when a coordinate and year appear as
    # both presence and absence.
    #
    # The presence-first sort is also the pipeline's row-order contract. GEE's libsvm
    # takes its positive class from the first label it sees in the training data, and
    # the SVM score flip in predict_gee_map() assumes that class is presence (1).
    # Emitting presences first here holds that invariant at the data boundary, so no
    # later stage has to re-sort. order() is stable, so the original order survives
    # within each class.
    rows_before <- nrow(result)

    key_cols <- setdiff(names(result), "present")
    result <- result[order(-result$present), ]
    result <- result[!duplicated(result[, key_cols, drop = FALSE]), ]
    rownames(result) <- NULL
    rows_after <- nrow(result)

    if (rows_before != rows_after) {
        sdm_warn(sprintf("%d row%s removed: duplicate or conflicting records (presence prioritized)",
                         rows_before - rows_after,
                         if ((rows_before - rows_after) == 1) "" else "s"))
    }

    cols_desc <- paste(names(result), collapse = ", ")
    prefix <- if (!is.null(label)) paste0(label, " data") else "Data"
    sdm_done(sprintf("%s ready: %d rows \u00d7 %d columns (%s)", prefix, nrow(result), ncol(result), cols_desc))

    return(result)
}
