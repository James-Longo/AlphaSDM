#' Install autoSDM Python Dependencies
#'
#' This function installs the required Python dependencies for autoSDM
#' using reticulate.
#'
#' @param method Installation method. Pass "auto" to let reticulate decide.
#' @param envname The name, or full path, of the environment in which Python packages are to be installed.
#' @export
install_autoSDM <- function(method = "virtualenv", envname = "r-autoSDM") {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("The 'reticulate' package is required to manage Python dependencies.")
    }

    packages <- c("earthengine-api")

    message("Installing Python dependencies for autoSDM into environment: ", envname)

    # Force creation if it doesn't exist to avoid reticulate's uv search
    if (!reticulate::virtualenv_exists(envname)) {
        message("Creating virtual environment...")
        reticulate::virtualenv_create(envname)
    }

    reticulate::py_install(packages, envname = envname, method = method, pip = TRUE)

    message("\nSuccess! GEE Python API installed.")
}

#' Ensure Google Earth Engine is Authenticated and Python is Setup (Internal)
#' @keywords internal
ensure_gee_authenticated <- function(project = NULL) {
    if (!requireNamespace("rgee", quietly = TRUE)) {
        stop("The 'rgee' package is required.")
    }

    # 0. Check session cache
    if (isTRUE(getOption("autoSDM.gee_initialized"))) {
        return(TRUE)
    }

    # 1. Automate Python Environment Selection
    # If the user hasn't set one, try to use our internal one
    if (!reticulate::py_available() || !reticulate::py_module_available("ee")) {
        py_path <- resolve_python_path()
        if (!is.null(py_path)) {
            try(reticulate::use_python(py_path, required = FALSE), silent = TRUE)
        }
    }

    # 2. Check if earthengine-api is actually installed
    if (!reticulate::py_module_available("ee")) {
        message("[autoSDM] GEE Python API not found. Attempting automatic installation...")
        install_autoSDM()
        # After install, we really should restart, but let's try to proceed
        py_path <- resolve_python_path()
        if (!is.null(py_path)) {
            try(reticulate::use_python(py_path, required = FALSE), silent = TRUE)
        }
    }

    config_dir <- file.path(Sys.getenv("HOME"), ".config", "autoSDM")
    config_file <- file.path(config_dir, "config.json")

    # 3. Determine Project ID
    # Prioritize the credentials file because it's what the user just authorized
    if (is.null(project) || project == "") {
        creds_file <- file.path(Sys.getenv("HOME"), ".config", "earthengine", "credentials")
        if (file.exists(creds_file)) {
            try({
                creds <- jsonlite::fromJSON(creds_file)
                project <- creds$project
            }, silent = TRUE)
        }
    }

    # Fallback 1: Environment variable
    if (is.null(project) || project == "") {
        project <- Sys.getenv("EARTHENGINE_PROJECT", unset = "")
    }

    # Fallback 2: Local autoSDM config
    if (is.null(project) || project == "") {
        if (file.exists(config_file)) {
            try({
                conf <- jsonlite::fromJSON(config_file)
                project <- conf$gee_project
            }, silent = TRUE)
        }
    }

    # 4. Try simple initialization
    message("Initializing Google Earth Engine...")
    success <- tryCatch(
        {
            if (!is.null(project) && project != "") {
                rgee::ee_Initialize(project = project, drive = FALSE, gcs = FALSE)
            } else {
                rgee::ee_Initialize(drive = FALSE, gcs = FALSE)
            }
            TRUE
        },
        error = function(e) {
            # If it fails, return the error object to check message
            return(e)
        }
    )

    if (isTRUE(success)) {
        options(autoSDM.gee_initialized = TRUE)
        return(TRUE)
    }

    # 5. Handle Failure
    err_msg <- success$message
    
    # Handle missing Python package error
    if (grepl("does not have the Python package", err_msg)) {
        message("\n[autoSDM] GEE Python dependencies are missing or misconfigured.")
        install_autoSDM()
        stop("autoSDM environment created. Please RESTART your R session to use the new environment.")
    }

    # Check for actual expiration vs project mismatch
    is_expired <- grepl("expired", err_msg, ignore.case = TRUE)
    is_project_error <- grepl("project", err_msg, ignore.case = TRUE) || grepl("permission", err_msg, ignore.case = TRUE)

    if (is_expired && !is_project_error) {
        message("\n[autoSDM] GEE credentials appear to be EXPIRED.")
        # Only offer cleanup if interactive
        if (interactive()) {
            confirm <- utils::askYesNo("Would you like to clear local credentials and re-authenticate?")
            if (isTRUE(confirm)) {
                 rgee::ee_clean_user_credentials()
                 rgee::ee_Authenticate()
                 message("\n[autoSDM] Credentials reset. PLEASE RESTART YOUR R SESSION for changes to take effect.")
                 stop("R session restart required after GEE credential reset.")
            }
        } else {
            stop("GEE credentials expired. Run rgee::ee_clean_user_credentials() and re-authenticate in an interactive session.")
        }
    } else if (is_project_error) {
        message("\n[autoSDM] GEE Initialization failed due to a Project ID or Permission issue.")
        message("Error: ", err_msg)
        if (is.null(project) || project == "") {
            message("\nTIP: Try specifying your Google Cloud Project ID explicitly:")
            message("autoSDM(..., gee_project = 'your-project-id')")
        }
        stop("GEE Project ID error. Please ensure your project is correct and GEE is enabled.")
    } else {
        # General failure - try one re-auth attempt
        message("\n[autoSDM] GEE initialization failed. Attempting setup...")
        rgee::ee_Authenticate()
        tryCatch({
            rgee::ee_Initialize(project = project, drive = FALSE, gcs = FALSE)
            options(autoSDM.gee_initialized = TRUE)
        }, error = function(e2) {
            stop("Final GEE initialization attempt failed. Please check your internet connection, Project ID, and GEE account status.\nError: ", e2$message)
        })
    }
    
    return(TRUE)
}

#' Hard Reset Google Earth Engine Credentials
#'
#' This function forcefully removes all local Google Earth Engine credentials
#' to resolve persistent authentication loops or "expired token" errors.
#' After running this, your next GEE-related call will prompt for a fresh browser authentication.
#'
#' @export
ee_hard_reset <- function() {
    if (!requireNamespace("rgee", quietly = TRUE)) stop("rgee package required.")
    
    message("Performing hard reset of Google Earth Engine credentials...")
    
    # Use rgee helper
    try(rgee::ee_clean_user_credentials(), silent = TRUE)
    
    # Manual target search
    ee_cfg <- file.path(Sys.getenv("HOME"), ".config", "earthengine")
    if (dir.exists(ee_cfg)) {
        message("Removing directory: ", ee_cfg)
        unlink(ee_cfg, recursive = TRUE)
    }
    
    message("Reset complete. Please restart your R session and run autoSDM again.")
}

#' Ensure dependencies are met (Internal)
#' @keywords internal
ensure_autoSDM_dependencies <- function() {
    if (!requireNamespace("rgee", quietly = TRUE)) {
        stop("The 'rgee' package is required.")
    }
    ensure_gee_authenticated()
    return(TRUE)
}

#' Resolve Python Path (Internal)
#' @keywords internal
resolve_python_path <- function(path = NULL) {
    if (!is.null(path)) return(path)
    if (reticulate::virtualenv_exists("r-autoSDM")) return(reticulate::virtualenv_python("r-autoSDM"))
    if (reticulate::py_available()) return(reticulate::py_config()$python)
    return(NULL)
}
