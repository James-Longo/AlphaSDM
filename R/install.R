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

#' Path to the live GEE credentials file (the one the ee client reads)
#' @keywords internal
.gee_live_cred_path <- function() {
    file.path(Sys.getenv("HOME"), ".config", "earthengine", "credentials")
}

#' Check if GEE credentials exist on disk
#' @keywords internal
.gee_credentials_exist <- function() {
    file.exists(.gee_live_cred_path())
}

#' Durable credential store (optional)
#'
#' On some platforms the home/config directory is ephemeral, wiped between
#' processes (sandboxes, some container and HPC-scratch setups), so the GEE
#' credentials the ee client writes to \code{~/.config/earthengine/credentials}
#' do not survive to the next run. If the environment variable
#' \code{ALPHASDM_GEE_CRED_STORE} points at a directory on a PERSISTENT
#' filesystem, AlphaSDM mirrors the credentials there and restores them at the
#' start of every session, so a one-time authentication persists.
#'
#' Returns the store directory (NULL if the feature is not enabled).
#' @keywords internal
.gee_cred_store <- function() {
    d <- Sys.getenv("ALPHASDM_GEE_CRED_STORE", unset = "")
    if (!nzchar(d)) return(NULL)
    d
}

#' Copy the live credentials into the durable store (after a successful auth)
#' @keywords internal
.gee_backup_credentials <- function() {
    store <- .gee_cred_store(); if (is.null(store)) return(invisible(FALSE))
    live <- .gee_live_cred_path(); if (!file.exists(live)) return(invisible(FALSE))
    dir.create(store, showWarnings = FALSE, recursive = TRUE)
    ok <- file.copy(live, file.path(store, "credentials"), overwrite = TRUE)
    proj_cfg <- file.path(Sys.getenv("HOME"), ".config", "AlphaSDM", "config.json")
    if (file.exists(proj_cfg)) file.copy(proj_cfg, file.path(store, "config.json"), overwrite = TRUE)
    invisible(ok)
}

#' Restore credentials from the durable store into the live location
#'
#' Runs at session start: if the live credentials are missing but a durable copy
#' exists, copy it into place (and the saved project id) so the ee client finds
#' a working token without re-authenticating. No-op when the store is unset or empty.
#' @keywords internal
.gee_restore_credentials <- function() {
    store <- .gee_cred_store(); if (is.null(store)) return(invisible(FALSE))
    src <- file.path(store, "credentials"); if (!file.exists(src)) return(invisible(FALSE))
    if (.gee_credentials_exist()) return(invisible(TRUE))   # live copy already present
    dir.create(dirname(.gee_live_cred_path()), showWarnings = FALSE, recursive = TRUE)
    ok <- file.copy(src, .gee_live_cred_path(), overwrite = TRUE)
    proj_src <- file.path(store, "config.json")
    if (file.exists(proj_src)) {
        dir.create(file.path(Sys.getenv("HOME"), ".config", "AlphaSDM"), showWarnings = FALSE, recursive = TRUE)
        file.copy(proj_src, file.path(Sys.getenv("HOME"), ".config", "AlphaSDM", "config.json"), overwrite = TRUE)
    }
    invisible(ok)
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

# ---- Internal Connection Helpers ----

#' Test whether a Python interpreter already has the earthengine-api module
#'
#' Cheap out-of-process probe: runs `python -c "find_spec('ee')"` and checks the
#' exit code. Does NOT import via reticulate (which would lock the session to
#' that interpreter), so it is safe to call against several candidates.
#' @keywords internal
.py_has_ee <- function(py) {
    if (is.null(py) || !nzchar(py) || !file.exists(py)) return(FALSE)
    code <- "import importlib.util, sys; sys.exit(0 if importlib.util.find_spec('ee') else 1)"
    ok <- tryCatch(
        suppressWarnings(system2(py, args = c("-c", shQuote(code)),
                                 stdout = FALSE, stderr = FALSE)),
        error = function(e) 1L)
    identical(as.integer(ok), 0L)
}

#' Find an existing Python that already provides earthengine-api
#'
#' Startup optimization: before building a fresh Python environment (slow, and
#' on some systems triggers a ~100 MB Miniconda download), look for an
#' interpreter that can already `import ee`. Checks, in order: the interpreter
#' reticulate is already bound to, EARTHENGINE_PYTHON / RETICULATE_PYTHON, the
#' active conda env, and reticulate's discovered default. Returns the path of
#' the first match, or NULL if none has earthengine-api.
#' @keywords internal
.find_ee_python <- function() {
    # Quiet wrapper: reticulate's config probes emit an alarming
    # "Unable to find conda binary" message on systems without conda even when
    # they still return a usable interpreter path. We only want the path.
    quiet <- function(expr) suppressWarnings(suppressMessages(tryCatch(
        utils::capture.output(val <- expr, type = "message"),
        error = function(e) NULL)))
    py_config_path <- function() {
        val <- NULL
        quiet({ cfg <- reticulate::py_config();          val <<- cfg$python })
        val
    }
    py_discover_path <- function() {
        val <- NULL
        quiet({ cfg <- reticulate::py_discover_config();  val <<- cfg$python })
        val
    }
    # Cheapest, most reliable candidates first; reticulate probes last.
    cands <- c(
        Sys.getenv("EARTHENGINE_PYTHON", ""),
        Sys.getenv("RETICULATE_PYTHON", ""),
        { cp <- Sys.getenv("CONDA_PREFIX", ""); if (nzchar(cp)) file.path(cp, "bin", "python") else "" },
        tryCatch(py_config_path(),   error = function(e) NULL),
        tryCatch(py_discover_path(), error = function(e) NULL)
    )
    cands <- unique(cands[!vapply(cands, is.null, logical(1))])
    cands <- cands[nzchar(cands)]
    for (p in cands) if (.py_has_ee(path.expand(p))) return(path.expand(p))
    NULL
}


#' Report the auth type stored in the persistent credentials file
#'
#' Returns "user account (OAuth)" for the normal personal-Google-account flow,
#' "service account" for a key file, or NA if no credentials are on disk.
#' @keywords internal
.gee_auth_type <- function() {
    f <- file.path(Sys.getenv("HOME"), ".config", "earthengine", "credentials")
    if (!file.exists(f)) return(NA_character_)
    info <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
    if (is.null(info)) return("unknown")
    if (!is.null(info$type) && identical(info$type, "service_account"))
        return("service account")
    if (!is.null(info$refresh_token)) return("user account (OAuth)")
    "unknown"
}

#' Detect a headless / no-browser session
#'
#' The zero-friction browser flow (auth_mode = "localhost") needs a web browser
#' on the same machine. That is always true on Windows and macOS. On Linux it
#' requires a display server (X11 or Wayland); a bare SSH shell, container, or
#' HPC node has neither, so localhost would hang. This lets setup_gee() pick a
#' working flow automatically instead of silently stalling.
#' @keywords internal
.gee_is_headless <- function() {
    sysname <- Sys.info()[["sysname"]]
    if (identical(sysname, "Windows") || identical(sysname, "Darwin")) return(FALSE)
    # Linux: viable only if a graphical display is present
    Sys.getenv("DISPLAY", "") == "" && Sys.getenv("WAYLAND_DISPLAY", "") == ""
}

#' Attempt a quiet GEE initialization, returning TRUE/FALSE
#'
#' Used as the idempotency probe: if existing on-disk credentials still work,
#' setup_gee() short-circuits and never prompts. Any failure (missing/expired
#' credentials, refresh error, bad project) returns FALSE so the caller can
#' fall through to interactive setup.
#' @keywords internal
.gee_try_init <- function(project) {
    .suppress_gee_deprecation_warnings()
    # rgee frequently re-signals real init failures as warnings, so during this
    # silent probe we treat any warning as a failure (fall through to setup).
    ok <- tryCatch({
        rgee::ee_Initialize(project = project, drive = FALSE, gcs = FALSE, quiet = TRUE)
        # Confirm the connection is actually live, not just locally configured:
        # a trivial server round-trip forces a token refresh + API call.
        ee <- reticulate::import("ee")
        identical(as.integer(ee$Number(1L)$getInfo()), 1L)
    }, error = function(e) FALSE, warning = function(w) FALSE)
    isTRUE(ok)
}

# ---- Exported Functions ----

#' Set Up Google Earth Engine for AlphaSDM (one-time)
#'
#' Connects AlphaSDM to Google Earth Engine using your personal Google account.
#' You only ever need to run this once per machine: it authenticates through a
#' single browser click and saves long-lived credentials, so every future R
#' session connects automatically with no further prompts.
#'
#' \strong{Before you start} you need a free Earth Engine account. Sign up at
#' \url{https://earthengine.google.com/signup/}. Earth Engine is free for
#' noncommercial, research, education, and nonprofit use. Registration links
#' your Google account to a Cloud project (its ID is what you pass as
#' \code{project}).
#'
#' \strong{The authentication is a browser click, not a code to paste.} On a
#' desktop or laptop, \code{setup_gee()} opens your browser, you click
#' \strong{Allow}, and the credential is captured automatically over a local
#' loopback port (\code{auth_mode = "localhost"}). Nothing is copied or pasted,
#' and the saved credentials do not expire with normal use.
#'
#' Re-running \code{setup_gee()} when you are already connected is a harmless
#' no-op; it detects the working credentials and returns immediately.
#'
#' @param project Google Cloud / Earth Engine project ID (e.g.,
#'   \code{"my-ee-project"}). If \code{NULL}, the saved project is reused, or you
#'   are prompted interactively.
#' @param force If \code{TRUE}, re-authenticate even if valid credentials
#'   already exist.
#' @param auth_mode Optional override for the Earth Engine authorization flow,
#'   passed straight to \code{ee$Authenticate()}. Leave \code{NULL} (recommended)
#'   to let AlphaSDM choose: \code{"localhost"} (one-click, no paste) on a
#'   machine with a browser, and \code{"notebook"} on a detected headless/remote
#'   session. Set \code{"notebook"} yourself to force the paste-a-code flow, or
#'   \code{"gcloud"} if you use the gcloud CLI.
#' @return Invisibly \code{TRUE} on success, or \code{FALSE} if a step (such as
#'   the Python install) requires you to restart R and re-run.
#' @export
setup_gee <- function(project = NULL, force = FALSE, auth_mode = NULL) {
    if (!requireNamespace("rgee", quietly = TRUE)) {
        stop("The 'rgee' package is required. Install it with: install.packages('rgee')")
    }

    # 1. Python environment.
    #    FAST PATH: if any interpreter reticulate can reach already has
    #    earthengine-api, bind to it and skip installation entirely. This avoids
    #    the slow rgee::ee_install() rebuild, and, on machines where reticulate
    #    cannot find a conda binary, the ~100 MB Miniconda download it falls back
    #    to. Most users already have a suitable Python (conda env, venv, or a
    #    system Python with `pip install earthengine-api`).
    ee_py <- .find_ee_python()
    if (!is.null(ee_py)) {
        Sys.setenv(EARTHENGINE_PYTHON = ee_py, RETICULATE_PYTHON = ee_py)
        sdm_done(sprintf("Using existing Python with earthengine-api: %s", ee_py))
    } else if (Sys.getenv("EARTHENGINE_PYTHON") != "" && !force) {
        # A Python was previously configured but earthengine-api isn't importable
        # from it: install just the package into it, no full env rebuild.
        sdm_section("Installing earthengine-api into your Python environment")
        ok <- tryCatch({ reticulate::py_install("earthengine-api", pip = TRUE); TRUE },
                       error = function(e) { sdm_warn(conditionMessage(e)); FALSE })
        if (!ok || !.py_has_ee(Sys.getenv("EARTHENGINE_PYTHON"))) {
            sdm_warn("Could not add earthengine-api to the current Python.")
            sdm_warn("Install it manually (pip install earthengine-api) or run setup_gee(force = TRUE).")
            return(invisible(FALSE))
        }
        sdm_done("earthengine-api installed")
    } else {
        # No suitable Python found: build rgee's environment as a last resort.
        # This is the only path that may download Miniconda; we surface that.
        sdm_section("Setting up a Python environment for Earth Engine")
        sdm_info("No Python with earthengine-api was found. Building one now")
        sdm_info("(first time only; this can take a few minutes).")
        rgee::ee_install(py_env = "rgee", confirm = FALSE)
        sdm_done("Python environment installed")
        sdm_warn("Please restart your R session, then run: AlphaSDM::setup_gee()")
        return(invisible(FALSE))
    }

    # Resolve project up front: argument > saved config > interactive prompt.
    if (is.null(project) || project == "") {
        project <- .read_saved_project()
    }

    # 2. Idempotency short-circuit: if we are NOT forcing and credentials on disk
    #    already produce a live connection, we are done; never prompt again.
    #    This is what makes setup a genuine one-time action.
    #    First restore from the durable store (if configured) in case the live
    #    ~/.config copy was wiped since last session (ephemeral-home platforms).
    .gee_restore_credentials()
    need_auth <- force || !.gee_credentials_exist()
    if (!force && .gee_credentials_exist()) {
        sdm_section("Checking existing Google Earth Engine connection")
        if (.gee_try_init(if (is.null(project) || project == "") NULL else trimws(project))) {
            options(AlphaSDM.gee_initialized = TRUE)
            .suppress_gee_deprecation_warnings()
            if (!is.null(project) && trimws(project) != "") .save_project(trimws(project))
            sdm_done(sprintf("Already connected to Earth Engine (%s). Nothing to do.",
                             .gee_auth_type()))
            return(invisible(TRUE))
        }
        # Credentials are present but do not produce a live connection (expired
        # or revoked). Force a fresh authentication rather than falling through
        # to ee_Initialize with stale credentials and a confusing error.
        sdm_info("Stored credentials are missing or expired; re-authenticating.")
        need_auth <- TRUE
    }

    # 3. Authenticate with the personal Google account.
    #    Default flow is auth_mode = "localhost": the browser opens, the user
    #    clicks Allow, and the token is captured automatically over a loopback
    #    port, with NO code to copy or paste, and the credentials are long-lived.
    #    On a detected headless/remote session (no browser reachable) we fall
    #    back to "notebook" so the flow still completes instead of hanging.
    if (need_auth) {
        chosen_mode <- auth_mode
        if (is.null(chosen_mode)) {
            if (.gee_is_headless()) {
                chosen_mode <- "notebook"
                sdm_section("Authenticating with Google Earth Engine (headless session)")
                sdm_info("No local browser detected. Open the printed URL on any device,")
                sdm_info("approve access, and paste the code back here.")
                sdm_info("Tip: on a desktop with a browser, this step is a single click, no paste.")
            } else {
                chosen_mode <- "localhost"
                sdm_section("Authenticating with Google Earth Engine")
                sdm_info("A browser window will open; click 'Allow'. No code to paste.")
            }
        } else {
            sdm_section(sprintf("Authenticating with Google Earth Engine (auth_mode = '%s')",
                                chosen_mode))
        }
        ee <- reticulate::import("ee")
        ee$Authenticate(auth_mode = chosen_mode, force = force)
        # Mirror the freshly written token to the durable store (if configured)
        # so it survives an ephemeral home directory on the next session.
        .gee_backup_credentials()
    }

    # 4. Get project ID (now that we are authenticated)
    if ((is.null(project) || project == "") && interactive()) {
        project <- readline("[AlphaSDM] Enter your Earth Engine / Cloud Project ID: ")
    }
    if (is.null(project) || trimws(project) == "") {
        stop("An Earth Engine / Google Cloud Project ID is required.\n",
             "It is created for you when you register at ",
             "https://earthengine.google.com/signup/ \n",
             "and is visible at https://console.cloud.google.com/", call. = FALSE)
    }
    project <- trimws(project)

    # 5. Verify the connection works
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

    # 6. Persist project ID for future sessions
    .save_project(project)
    .gee_backup_credentials()   # mirror creds + project id to the durable store, if configured
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

#' Report the Google Earth Engine Connection Status
#'
#' Prints a quick diagnostic of how AlphaSDM is connected to Earth Engine:
#' which Python environment is bound, whether saved credentials exist and are
#' the personal-account (OAuth) type, which project is configured, and whether a
#' live connection succeeds. Use it to confirm setup or to troubleshoot.
#'
#' Note: this reports the *connection*. To monitor running server-side export
#' tasks, use \code{\link{sdm_gee_status}} instead.
#'
#' @param check_live If \code{TRUE} (default), perform a small server round-trip
#'   to confirm the credentials actually work, not just that they are on disk.
#' @return Invisibly, a named list of the status fields.
#' @export
gee_status <- function(check_live = TRUE) {
    sdm_section("AlphaSDM: Google Earth Engine connection")

    py_env  <- Sys.getenv("EARTHENGINE_PYTHON", unset = "")
    creds   <- .gee_credentials_exist()
    atype   <- .gee_auth_type()
    project <- .read_saved_project()
    if (is.null(project) || project == "")
        project <- Sys.getenv("EARTHENGINE_PROJECT", unset = "")

    mark <- function(ok) if (isTRUE(ok)) "OK  " else "MISSING"

    sdm_info(sprintf("[%s] Python env : %s", mark(py_env != ""),
                     if (py_env == "") "not set (run setup_gee once)" else py_env), indent = 1L)
    sdm_info(sprintf("[%s] Credentials: %s", mark(creds),
                     if (creds) sprintf("present, %s", atype) else "none on disk"), indent = 1L)
    sdm_info(sprintf("[%s] Project    : %s", mark(project != ""),
                     if (project == "") "not set" else project), indent = 1L)

    live <- NA
    if (check_live && creds) {
        live <- .gee_try_init(if (project == "") NULL else project)
        sdm_info(sprintf("[%s] Live check : %s", mark(live),
                         if (isTRUE(live)) "connected" else "could not reach Earth Engine"),
                 indent = 1L)
    }

    if (!creds) {
        sdm_info("Not connected yet. Run: AlphaSDM::setup_gee(project = 'your-project-id')")
    } else if (isTRUE(live) || !check_live) {
        sdm_done("Earth Engine is set up. No action needed.")
    } else {
        sdm_warn("Credentials exist but the live check failed. Try: setup_gee(force = TRUE)")
    }

    invisible(list(python_env = py_env, credentials = creds, auth_type = atype,
                   project = project, live = live))
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

    # Restore credentials from the durable store if the live copy was wiped
    # (ephemeral-home platforms); a no-op when the store is unset or the live
    # copy is already present.
    .gee_restore_credentials()

    # Resolve project: argument > saved config > env var
    if (is.null(project) || project == "") {
        project <- .read_saved_project()
    }
    if (is.null(project) || project == "") {
        project <- Sys.getenv("EARTHENGINE_PROJECT", unset = "")
    }
    if (project == "") project <- NULL

    # Make sure reticulate is pointed at a Python that actually has
    # earthengine-api. reticulate's default discovery can bind to an unrelated
    # interpreter (e.g. a project-local .venv) that lacks `ee`; setting the env
    # var before the first import avoids that. We set the env var rather than
    # calling use_python() to sidestep reticulate's conda-binary lookup.
    if (Sys.getenv("EARTHENGINE_PYTHON") == "") {
        ee_py <- .find_ee_python()
        if (!is.null(ee_py)) Sys.setenv(EARTHENGINE_PYTHON = ee_py, RETICULATE_PYTHON = ee_py)
    }

    # Suppress GEE Python DeprecationWarnings BEFORE initializing: the GEE
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

    # If it failed, stop and tell the user to run the one-time setup
    stop("Not connected to Google Earth Engine.\n",
         "Run the one-time setup (a single browser click, no code to paste):\n",
         "    AlphaSDM::setup_gee(project = 'your-project-id')\n",
         "Check status any time with AlphaSDM::gee_status().\n",
         "Original error: ", result$message)
}
