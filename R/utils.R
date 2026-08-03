# ============================================================
#  AlphaSDM: Internal Logging & Progress Utilities (base R)
# ============================================================

# ---- Package-level state ----
.alphasdm_env <- new.env(parent = emptyenv())
.alphasdm_env$verbose <- TRUE
.alphasdm_env$standardization_active <- FALSE
.alphasdm_env$standardization_info_printed <- FALSE

# ---- Verbosity control ----

#' Control AlphaSDM console output verbosity
#'
#' @param verbose Logical. Set FALSE to suppress all non-error messages.
#' @export
sdm_verbose <- function(verbose = TRUE) {
  .alphasdm_env$verbose <- isTRUE(verbose)
  invisible(verbose)
}

# ---- Internal helper ----

.sdm_msg <- function(text) {
  if (!isTRUE(.alphasdm_env$verbose)) return(invisible(NULL))
  message(text)
}

# ============================================================
#  Public logging API
# ============================================================

#' Print a top-level section header
#' @param title Section title string.
#' @noRd
sdm_section <- function(title) {
  rule <- strrep("\u2500", max(0, 60 - nchar(title) - 2))
  .sdm_msg(sprintf("\n\u250c\u2500 %s %s", title, rule))
}

#' Print an informational status line
#' @param msg Message string.
#' @param indent Integer indent level (0 = top-level, 1 = nested).
#' @noRd
sdm_info <- function(msg, indent = 0L) {
  pad <- strrep("  ", indent)
  .sdm_msg(sprintf("%s  \u2023 %s", pad, msg))
}

#' Print a success / completion note
#' @param msg Message string.
#' @param indent Integer indent level.
#' @noRd
sdm_done <- function(msg, indent = 0L) {
  pad <- strrep("  ", indent)
  .sdm_msg(sprintf("%s  \u2714 %s", pad, msg))
}

#' Print a warning / advisory note (non-fatal)
#' @param msg Message string.
#' @param indent Integer indent level.
#' @noRd
sdm_warn <- function(msg, indent = 0L) {
  pad <- strrep("  ", indent)
  .sdm_msg(sprintf("%s  ! %s", pad, msg))
}

#' Start a named progress tracker and return a handle
#' @param name Label for this operation.
#' @return A list handle used with \code{sdm_progress_done}.
#' @noRd
sdm_progress_start <- function(name) {
  list(name = name, start = proc.time()[["elapsed"]])
}

#' Finish a progress tracker and print elapsed time
#' @param handle Handle returned by \code{sdm_progress_start}.
#' @noRd
sdm_progress_done <- function(handle) {
  if (is.null(handle) || !isTRUE(.alphasdm_env$verbose)) return(invisible(NULL))
  elapsed <- proc.time()[["elapsed"]] - handle$start
  .sdm_msg(sprintf("  \u2714 %s complete [%.1fs]", handle$name, elapsed))
  invisible(NULL)
}

#' Reset the per-run console state
#'
#' `format_data()` prints its section header only once per run. Registering this
#' with `on.exit()` means the state is cleared even when a run aborts, which would
#' otherwise leave the next `format_data()` in the session silently header-less.
#' @noRd
reset_sdm_run_state <- function() {
  .alphasdm_env$standardization_active <- FALSE
  .alphasdm_env$standardization_info_printed <- FALSE
  invisible(NULL)
}

#' Close out a run: elapsed-time line under a closing section header
#' @noRd
sdm_finish <- function(t_start, title) {
  sdm_section(title)
  sdm_done(sprintf("Total elapsed time [%.1fs]", proc.time()[["elapsed"]] - t_start))
  .sdm_msg("")
}
