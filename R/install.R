# ---- Python and the saved project ----

.onLoad <- function(libname, pkgname) {
    # Declare the Python dependency. reticulate provisions it the first time
    # Python starts, unless the user has chosen an interpreter of their own.
    reticulate::py_require("earthengine-api")
}

#' Path to the file holding the saved Earth Engine project id
#' @noRd
.alphasdm_config_file <- function() {
    file.path(tools::R_user_dir("AlphaSDM", "config"), "config.json")
}

#' Read the saved Earth Engine project id
#'
#' Falls back to the file used before version 0.2.0, which lived under
#' ~/.config/AlphaSDM; it is only read, never written.
#' @noRd
.read_saved_project <- function() {
    files <- c(.alphasdm_config_file(),
               file.path(Sys.getenv("HOME"), ".config", "AlphaSDM", "config.json"))
    for (f in files[file.exists(files)]) {
        p <- tryCatch(jsonlite::fromJSON(f)$gee_project, error = function(e) NULL)
        if (!is.null(p) && nzchar(p)) return(p)
    }
    NULL
}

#' Save the Earth Engine project id for later sessions
#' @noRd
.save_project <- function(project) {
    f <- .alphasdm_config_file()
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    jsonlite::write_json(list(gee_project = project), f, auto_unbox = TRUE)
}

#' Resolve the Earth Engine project: the argument, then the saved project, then
#' the EARTHENGINE_PROJECT environment variable. NULL when none is set.
#' @noRd
.resolve_project <- function(project = NULL) {
    if (is.null(project) || !nzchar(trimws(project))) project <- .read_saved_project()
    if (is.null(project) || !nzchar(project)) project <- Sys.getenv("EARTHENGINE_PROJECT")
    if (nzchar(project)) trimws(project) else NULL
}

#' Path of the credentials file the Earth Engine client writes
#' @noRd
.gee_cred_path <- function() {
    tryCatch(reticulate::import("ee")$oauth$get_credentials_path(),
             error = function(e) file.path(Sys.getenv("HOME"), ".config",
                                           "earthengine", "credentials"))
}

#' Report the kind of credentials on disk: "user account (OAuth)", "service
#' account", "unknown", or NA when there are none.
#' @noRd
.gee_auth_type <- function() {
    f <- .gee_cred_path()
    if (!file.exists(f)) return(NA_character_)
    info <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
    if (is.null(info)) return("unknown")
    if (identical(info$type, "service_account")) return("service account")
    if (!is.null(info$refresh_token)) return("user account (OAuth)")
    "unknown"
}

#' Whether this session has no local browser for the one-click sign-in
#'
#' Windows and macOS always have one; Linux needs a display server, which a bare
#' SSH shell, container or HPC node lacks.
#' @noRd
.gee_is_headless <- function() {
    if (Sys.info()[["sysname"]] %in% c("Windows", "Darwin")) return(FALSE)
    !nzchar(Sys.getenv("DISPLAY")) && !nzchar(Sys.getenv("WAYLAND_DISPLAY"))
}

#' Silence Earth Engine's notices about deprecated catalog assets
#'
#' The client re-enables these warnings during Initialize(), so the only
#' reliable fix is to replace the function that issues them. AlphaSDM uses no
#' deprecated assets.
#' @noRd
.suppress_gee_deprecation_warnings <- function() {
    try(reticulate::py_run_string(paste0(
        "try:\n",
        "    import ee.deprecation as _ee_dep\n",
        "    _ee_dep._IssueAssetDeprecationWarning = lambda asset: None\n",
        "    _ee_dep.InitializeDeprecatedAssets = lambda: None\n",
        "except Exception:\n",
        "    pass")), silent = TRUE)
}

#' Initialise Earth Engine for a project, with errors a user can act on
#' @noRd
.gee_init <- function(project) {
    ee <- tryCatch(reticulate::import("ee"), error = function(e)
        stop("The Earth Engine Python client (earthengine-api) is not available.\n",
             "If you set RETICULATE_PYTHON to your own Python, install it there with\n",
             "  pip install earthengine-api\n",
             "Original error: ", conditionMessage(e), call. = FALSE))
    .suppress_gee_deprecation_warnings()
    tryCatch(ee$Initialize(project = project), error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("Billing is disabled", msg, ignore.case = TRUE))
            stop(sprintf(paste0("Billing is not enabled for project '%s'. Even free Earth ",
                                "Engine use needs a linked billing account: ",
                                "https://console.cloud.google.com/billing/linkedaccount?project=%s"),
                         project, project), call. = FALSE)
        if (grepl("not registered|not found", msg, ignore.case = TRUE))
            stop(sprintf(paste0("Project '%s' is not registered for Earth Engine. Register ",
                                "at https://code.earthengine.google.com/register"), project),
                 call. = FALSE)
        stop(msg, call. = FALSE)
    })
    .suppress_gee_deprecation_warnings()
    .alphasdm_env$gee_initialized <- TRUE
    invisible(TRUE)
}

#' Initialise and make one round trip, returning TRUE or FALSE
#'
#' Confirms that credentials actually work, not only that they are on disk.
#' @noRd
.gee_try_init <- function(project) {
    isTRUE(tryCatch({
        .gee_init(project)
        identical(as.integer(reticulate::import("ee")$Number(1L)$getInfo()), 1L)
    }, error = function(e) FALSE))
}

# ---- Exported functions ----

#' Connect AlphaSDM to Google Earth Engine (one-time)
#'
#' Signs in to Google Earth Engine with your own Google account and saves your
#' project ID, so later R sessions connect on their own. Run it once per
#' machine; running it again when already connected does nothing.
#'
#' You need a free Earth Engine account first: register at
#' \url{https://earthengine.google.com/signup/}. Earth Engine is free for
#' noncommercial, research, education and nonprofit use, and registration gives
#' you the Cloud project ID to pass as \code{project}.
#'
#' Signing in opens your browser, where you click \strong{Allow}. The Earth
#' Engine client stores the resulting credentials in its own configuration
#' folder, as it does for every tool that uses Earth Engine. AlphaSDM saves only
#' the project ID, in \code{tools::R_user_dir("AlphaSDM", "config")}.
#'
#' The Earth Engine Python client (\code{earthengine-api}) is provided through
#' \pkg{reticulate}, which sets up a Python environment for it the first time
#' it is needed, unless you have pointed reticulate at a Python of your own.
#'
#' @param project Earth Engine (Google Cloud) project ID, for example
#'   \code{"my-ee-project"}. If \code{NULL}, the saved project is used, or you
#'   are asked for one.
#' @param force If \code{TRUE}, sign in again even if valid credentials exist.
#' @param auth_mode Earth Engine sign-in flow, passed to
#'   \code{ee.Authenticate()}. Leave \code{NULL} to use \code{"localhost"} (a
#'   browser click) where a browser is available and \code{"notebook"} (paste a
#'   code) where it is not. \code{"gcloud"} uses the gcloud command-line tool.
#' @return Invisibly, \code{TRUE} once connected.
#' @examples
#' \dontrun{
#' # Needs an Earth Engine account and an interactive session.
#' setup_gee(project = "my-ee-project")
#' }
#' @export
setup_gee <- function(project = NULL, force = FALSE, auth_mode = NULL) {
    project <- .resolve_project(project)

    # Already signed in: confirm the connection works and stop there.
    if (!force && file.exists(.gee_cred_path()) && !is.null(project)) {
        sdm_section("Checking the existing Earth Engine connection")
        if (.gee_try_init(project)) {
            .save_project(project)
            sdm_done(sprintf("Already connected to Earth Engine (%s). Nothing to do.",
                             .gee_auth_type()))
            return(invisible(TRUE))
        }
        sdm_info("The stored credentials no longer work; signing in again.")
    }

    if (!interactive())
        stop("Signing in to Earth Engine needs an interactive R session. Run ",
             "setup_gee() once in an interactive session on this machine.", call. = FALSE)

    if (is.null(auth_mode))
        auth_mode <- if (.gee_is_headless()) "notebook" else "localhost"
    sdm_section("Signing in to Google Earth Engine")
    if (identical(auth_mode, "notebook")) {
        sdm_info("Open the printed URL on any device, approve access, and paste the code back here.")
    } else if (identical(auth_mode, "localhost")) {
        sdm_info("A browser window will open; click 'Allow'.")
    }
    reticulate::import("ee")$Authenticate(auth_mode = auth_mode, force = TRUE)

    if (is.null(project))
        project <- trimws(readline("Earth Engine (Google Cloud) project ID: "))
    if (!nzchar(project))
        stop("An Earth Engine project ID is required. Registering at ",
             "https://earthengine.google.com/signup/ creates one; it is listed at ",
             "https://console.cloud.google.com/", call. = FALSE)

    sdm_section("Verifying the Earth Engine connection")
    .gee_init(project)
    .save_project(project)
    sdm_done(sprintf("Connected. Project '%s' saved for future sessions.", project))
    invisible(TRUE)
}

#' Forget the saved Earth Engine project, and optionally sign out
#'
#' Removes the project ID AlphaSDM saved. In an interactive session it then
#' offers to delete the Earth Engine sign-in credentials as well; those are
#' shared by every tool on this computer that uses Earth Engine, so they are
#' kept unless you agree. Run \code{\link{setup_gee}} afterwards to reconnect.
#'
#' @return Invisibly, \code{TRUE}.
#' @examples
#' \dontrun{
#' clear_gee_credentials()
#' }
#' @export
clear_gee_credentials <- function() {
    config_file <- .alphasdm_config_file()
    if (file.exists(config_file)) {
        unlink(config_file)
        sdm_done(sprintf("Removed the saved project: %s", config_file))
    }
    cred <- .gee_cred_path()
    if (file.exists(cred) && interactive() && isTRUE(utils::askYesNo(paste0(
        "Also sign out of Earth Engine on this computer? This deletes ", cred,
        ", which other Earth Engine tools use too."), default = FALSE))) {
        unlink(cred)
        sdm_done(sprintf("Removed: %s", cred))
    }
    .alphasdm_env$gee_initialized <- NULL
    invisible(TRUE)
}

#' Report the Google Earth Engine connection status
#'
#' Prints whether the Earth Engine client is available, whether sign-in
#' credentials exist and of which kind, which project is configured, and
#' whether a live connection succeeds. To monitor running Earth Engine tasks,
#' use \code{\link{gee_tasks}} instead.
#'
#' @param check_live If \code{TRUE} (default), make a small request to confirm
#'   that the credentials work, not only that they are on disk.
#' @return Invisibly, a named list of the status fields.
#' @examples
#' \dontrun{
#' gee_status()
#' }
#' @export
gee_status <- function(check_live = TRUE) {
    sdm_section("AlphaSDM: Google Earth Engine connection")
    client  <- tryCatch(reticulate::import("ee")$`__version__`, error = function(e) NA_character_)
    creds   <- file.exists(.gee_cred_path())
    atype   <- .gee_auth_type()
    project <- .resolve_project()
    mark    <- function(ok) if (isTRUE(ok)) "OK  " else "MISSING"

    sdm_info(sprintf("[%s] Python client: %s", mark(!is.na(client)),
                     if (is.na(client)) "earthengine-api not available"
                     else paste("earthengine-api", client)), indent = 1L)
    sdm_info(sprintf("[%s] Credentials  : %s", mark(creds),
                     if (creds) atype else "none; run setup_gee()"), indent = 1L)
    sdm_info(sprintf("[%s] Project      : %s", mark(!is.null(project)),
                     if (is.null(project)) "not set" else project), indent = 1L)

    live <- NA
    if (check_live && creds && !is.null(project)) {
        live <- .gee_try_init(project)
        sdm_info(sprintf("[%s] Live check   : %s", mark(live),
                         if (live) "connected" else "could not reach Earth Engine"), indent = 1L)
    }
    if (!creds || is.null(project)) {
        sdm_info("Not connected yet. Run: setup_gee(project = 'your-project-id')")
    } else if (isTRUE(live) || !check_live) {
        sdm_done("Earth Engine is set up.")
    } else {
        sdm_warn("Credentials exist but the live check failed. Try setup_gee(force = TRUE).")
    }
    invisible(list(client = client, credentials = creds, auth_type = atype,
                   project = project, live = live))
}

# ---- Authentication gate ----

#' Make sure Earth Engine is initialised (internal)
#'
#' Called at the top of every user-facing function that talks to Earth Engine.
#' @param project Optional project ID override.
#' @noRd
ensure_gee_authenticated <- function(project = NULL) {
    if (isTRUE(.alphasdm_env$gee_initialized)) return(invisible(TRUE))
    project <- .resolve_project(project)
    res <- tryCatch(.gee_init(project), error = function(e) e)
    if (inherits(res, "error"))
        stop("Not connected to Google Earth Engine. Run the one-time setup:\n",
             "    setup_gee(project = 'your-project-id')\n",
             "and check it with gee_status().\nOriginal error: ",
             conditionMessage(res), call. = FALSE)
    invisible(TRUE)
}
