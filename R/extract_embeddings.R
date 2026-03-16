#' Extract Alpha Earth Embeddings
#'
#' This function extracts the 64-dimensional Alpha Earth satellite embeddings for a set of coordinates.
#' Standardizes extraction using Google Earth Engine.
#'
#' @param df A data frame containing at least `longitude`, `latitude`, and `year`.
#' @param scale Optional. Resolution in meters for extraction. Defaults to 10.
#' @param python_path Optional. Path to Python executable.
#' @param gee_project Optional. Google Cloud Project ID. If NULL, uses local config or environment variables.
#' @return A data frame containing the original columns plus 64 embedding columns (A00-A63).
#' @export
extract_embeddings <- function(df, scale = 10, python_path = NULL, gee_project = NULL) {
    # 1. Config & Validation
    python_path <- resolve_python_path(python_path)
    if (is.null(python_path)) stop("Python path not found. Please run install_autoSDM().")

    # 2. GEE Auth & Project Selection (Option B support)
    ensure_gee_authenticated(project = gee_project)

    # Fallback to local config if still NULL
    if (is.null(gee_project)) {
        config_file <- file.path(Sys.getenv("HOME"), ".config", "autoSDM", "config.json")
        if (file.exists(config_file)) {
            try(
                {
                    conf <- jsonlite::fromJSON(config_file)
                    gee_project <- conf$gee_project
                },
                silent = TRUE
            )
        }
    }

    # 3. Execution
    tmp_in <- tempfile(fileext = ".csv")
    tmp_out <- tempfile(fileext = ".csv")
    write.csv(df, tmp_in, row.names = FALSE)

    args <- c(
        "-m", "autoSDM.cli", "background",
        "--input", shQuote(tmp_in),
        "--output", shQuote(tmp_out),
        "--scale", scale,
        "--count", 0 # Mode used for extraction only
    )

    if (!is.null(gee_project)) args <- c(args, "--project", shQuote(gee_project))

    message(sprintf("Extracting Alpha Earth embeddings for %d points...", nrow(df)))
    status <- system2(python_path, args = args, stdout = "", stderr = "")

    if (status != 0) {
        unlink(tmp_in)
        unlink(tmp_out)
        stop("Embedding extraction failed.")
    }

    # 4. Load Results
    res <- read.csv(tmp_out)
    unlink(tmp_in)
    unlink(tmp_out)

    return(res)
}

#' Initialize Earth Engine with Service Account
#'
#' @param json_path Path to the service account JSON key. Defaults to GEE_SERVICE_ACCOUNT_KEY env var.
#' @param venv_path Path to the Python virtual environment.
#' @export
#' @keywords internal
ee_auth_service <- function(json_path = Sys.getenv("GEE_SERVICE_ACCOUNT_KEY"), venv_path = NULL) {
    if (json_path == "" || !file.exists(json_path)) {
        stop("Valid service account JSON key path required.")
    }

    # We can't easily switch projects once ee is initialized in reticulate without restarting,
    # but we can try to rely on the CLI's --key argument for subsequent calls.
    # This function primarily exists to register the key path for the session.
    assign("sa_json_key", json_path, envir = .GlobalEnv)
    message("Service account registered for this session.")
}
