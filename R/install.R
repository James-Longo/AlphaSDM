# ---- Internal Helpers ----

#' Read saved GEE project ID from AlphaSDM config
#' @keywords internal
.read_saved_project <- function() {
    f <- file.path(Sys.getenv("HOME"), ".config", "AlphaSDM", "config.json")
    if (!file.exists(f)) return(NULL)
    tryCatch(jsonlite::fromJSON(f)$gee_project, error = function(e) NULL)
}

#' Save GEE project ID to AlphaSDM config
#' @keywords internal
.save_project <- function(project) {
    dir <- file.path(Sys.getenv("HOME"), ".config", "AlphaSDM")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    jsonlite::write_json(list(gee_project = project),
                         file.path(dir, "config.json"), auto_unbox = TRUE)
}

#' Check if GEE credentials exist on disk
#' @keywords internal
.gee_credentials_exist <- function() {
    file.exists(file.path(Sys.getenv("HOME"), ".config", "earthengine", "credentials"))
}

#' Suppress Python DeprecationWarnings from the GEE client library
#'
#' GEE's deprecation.py actively defeats warnings.filterwarnings() by calling
#' _UnfilterDeprecationWarnings() during ee.Initialize(). The only reliable fix
#' is to monkey-patch _IssueAssetDeprecationWarning to a no-op before Initialize
#' runs. These are server-side catalog notices about deprecated GEE assets that
#' are entirely irrelevant to AlphaSDM's embedding pipeline.
#' @keywords internal
.suppress_gee_deprecation_warnings <- function() {
    tryCatch(
        reticulate::py_run_string(paste0(
            "try:\n",
            "    import ee.deprecation as _ee_dep\n",
            "    _ee_dep._IssueAssetDeprecationWarning = lambda asset: None\n",
            "    _ee_dep.InitializeDeprecatedAssets = lambda: None\n",
            "except Exception:\n",
            "    pass"
        )),
        error = function(e) NULL
    )
}

# ---- Exported Functions ----

#' Set Up Google Earth Engine for AlphaSDM
#'
#' One-time interactive setup. Installs Python dependencies via rgee,
#' authenticates with GEE through a browser OAuth flow, and saves your
#' project ID locally. After running this once, all AlphaSDM functions
#' connect to GEE automatically without further prompts.
#'
#' @param project Google Cloud Project ID (e.g., "my-ee-project").
#'   If NULL, you will be prompted to enter it interactively.
#' @param force If TRUE, re-run the full setup even if already configured.
#' @export
setup_gee <- function(project = NULL, force = FALSE) {
    if (!requireNamespace("rgee", quietly = TRUE)) {
        stop("The 'rgee' package is required. Install it with: install.packages('rgee')")
    }

    # 1. Python environment: use rgee's own installer
    #    Creates a virtualenv, installs earthengine-api + numpy,
    #    and writes EARTHENGINE_PYTHON / EARTHENGINE_ENV to ~/.Renviron
    if (force || Sys.getenv("EARTHENGINE_PYTHON") == "") {
        sdm_section("Installing Python environment")
        rgee::ee_install(py_env = "rgee", confirm = FALSE)
        sdm_done("Python environment installed")
        sdm_warn("Please restart your R session, then run: AlphaSDM::setup_gee()")
        return(invisible(FALSE))
    }

    # 2. Authenticate (browser OAuth - creates persistent token)
    #    Use Python API directly with auth_mode="localhost" for zero-friction:
    #    browser opens, user clicks Allow, token auto-captured via local redirect.
    if (force || !.gee_credentials_exist()) {
        sdm_section("Authenticating with Google Earth Engine")
        sdm_info("A browser window will open — click 'Allow' to grant access")
        ee <- reticulate::import("ee")
        ee$Authenticate(auth_mode = "localhost")
    }

    # 3. Get project ID
    if (is.null(project) || project == "") {
        project <- .read_saved_project()
    }
    if ((is.null(project) || project == "") && interactive()) {
        project <- readline("[AlphaSDM] Enter your Google Cloud Project ID: ")
    }
    if (is.null(project) || trimws(project) == "") {
        stop("A Google Cloud Project ID is required. ",
             "Get one at https://console.cloud.google.com/ ",
             "and enable the Earth Engine API.")
    }
    project <- trimws(project)

    # 4. Verify the connection works
    sdm_section("Verifying GEE connection")
    tryCatch({
        rgee::ee_Initialize(project = project, drive = FALSE, gcs = FALSE, quiet = TRUE)
    }, warning = function(w) {
        # rgee wraps the real GEE error as a warning; catch and re-raise clearly
        msg <- conditionMessage(w)
        if (grepl("Billing is disabled", msg, ignore.case = TRUE)) {
            stop(sprintf(
                "Billing is not enabled for project '%s'.\n",
                project,
                "Even for free-tier GEE usage, a billing account must be linked.\n",
                "Enable it at: https://console.cloud.google.com/billing/linkedaccount?project=",
                project
            ), call. = FALSE)
        }
        if (grepl("not registered", msg, ignore.case = TRUE) || grepl("not found", msg, ignore.case = TRUE)) {
            stop(sprintf(
                "Project '%s' is not registered for Earth Engine.\n%s%s",
                project,
                "Register at: https://code.earthengine.google.com/register\n",
                "Make sure the Earth Engine API is enabled in your Google Cloud Console."
            ), call. = FALSE)
        }
        # For other warnings, re-raise the original error
        warning(w)
    }, error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("Billing is disabled", msg, ignore.case = TRUE)) {
            stop(sprintf(
                "Billing is not enabled for project '%s'.\n%s%s%s",
                project,
                "Even for free-tier GEE usage, a billing account must be linked.\n",
                "Enable it at: https://console.cloud.google.com/billing/linkedaccount?project=",
                project
            ), call. = FALSE)
        }
        stop("GEE initialization failed for project '", project, "'.\nError: ", msg, call. = FALSE)
    })

    # 5. Persist project ID for future sessions
    .save_project(project)
    options(AlphaSDM.gee_initialized = TRUE)
    .suppress_gee_deprecation_warnings()

    sdm_done(sprintf("Setup complete! Project '%s' saved for future sessions.", project))
    invisible(TRUE)
}

#' Clear All GEE Credentials and Configuration
#'
#' Removes all locally stored GEE credentials and the saved project ID.
#' After calling this, you will need to run \code{\link{setup_gee}} again.
#'
#' @export
clear_gee_credentials <- function() {
    # Clear rgee's internal credentials
    if (requireNamespace("rgee", quietly = TRUE)) {
        try(rgee::ee_clean_user_credentials(), silent = TRUE)
    }

    # Clear earthengine config directory
    ee_cfg <- file.path(Sys.getenv("HOME"), ".config", "earthengine")
    if (dir.exists(ee_cfg)) {
        unlink(ee_cfg, recursive = TRUE)
        sdm_done(sprintf("Removed: %s", ee_cfg))
    }

    # Clear our saved config
    config_file <- file.path(Sys.getenv("HOME"), ".config", "AlphaSDM", "config.json")
    if (file.exists(config_file)) {
        unlink(config_file)
        sdm_done(sprintf("Removed: %s", config_file))
    }

    # Clear session cache
    options(AlphaSDM.gee_initialized = NULL)

    sdm_done("All GEE credentials cleared. Run setup_gee() to reconfigure.")
    invisible(TRUE)
}

# ---- Internal Authentication Gate ----

#' Ensure GEE is authenticated and initialized (internal)
#'
#' Called at the top of every user-facing function. Silently connects
#' if credentials exist, or triggers setup_gee() if interactive.
#'
#' @param project Optional project ID override.
#' @keywords internal
ensure_gee_authenticated <- function(project = NULL) {
    # Session cache - already initialized this R session
    if (isTRUE(getOption("AlphaSDM.gee_initialized"))) {
        return(TRUE)
    }

    if (!requireNamespace("rgee", quietly = TRUE)) {
        stop("The 'rgee' package is required. Install it with: install.packages('rgee')")
    }

    # Resolve project: argument > saved config > env var
    if (is.null(project) || project == "") {
        project <- .read_saved_project()
    }
    if (is.null(project) || project == "") {
        project <- Sys.getenv("EARTHENGINE_PROJECT", unset = "")
    }
    if (project == "") project <- NULL

    # Suppress GEE Python DeprecationWarnings BEFORE initializing — the GEE
    # client fires these during ee_Initialize itself (server-side catalog audit).
    # EARTHENGINE_PYTHON in .Renviron ensures rgee uses the correct virtualenv;
    # we do not override use_python() here to avoid re-initialization conflicts.
    .suppress_gee_deprecation_warnings()

    # Try to initialize (rgee handles expired token retry internally)
    result <- tryCatch({
        rgee::ee_Initialize(project = project, drive = FALSE, gcs = FALSE, quiet = TRUE)
        TRUE
    }, error = function(e) e)

    if (isTRUE(result)) {
        options(AlphaSDM.gee_initialized = TRUE)
        # Re-apply AFTER Initialize: GEE's deprecation.py calls _UnfilterDeprecationWarnings()
        # internally which inserts a 'default' filter for 'ee.deprecation', overriding our
        # pre-init filter. Re-applying here prepends a new 'ignore' that takes precedence.
        .suppress_gee_deprecation_warnings()
        return(TRUE)
    }

    # If it failed, stop and tell the user to run setup
    stop("GEE not configured or initialization failed.\n",
         "Please run AlphaSDM::setup_gee(project = 'your-project-id') in an interactive session first.\n",
         "Original error: ", result$message)
}
