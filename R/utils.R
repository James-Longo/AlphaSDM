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
#' @keywords internal
sdm_section <- function(title) {
  rule <- strrep("\u2500", max(0, 60 - nchar(title) - 2))
  .sdm_msg(sprintf("\n\u250c\u2500 %s %s", title, rule))
}

#' Print an informational status line
#' @param msg Message string.
#' @param indent Integer indent level (0 = top-level, 1 = nested).
#' @keywords internal
sdm_info <- function(msg, indent = 0L) {
  pad <- strrep("  ", indent)
  .sdm_msg(sprintf("%s  \u2023 %s", pad, msg))
}

#' Print a success / completion note
#' @param msg Message string.
#' @param indent Integer indent level.
#' @keywords internal
sdm_done <- function(msg, indent = 0L) {
  pad <- strrep("  ", indent)
  .sdm_msg(sprintf("%s  \u2714 %s", pad, msg))
}

#' Print a warning / advisory note (non-fatal)
#' @param msg Message string.
#' @param indent Integer indent level.
#' @keywords internal
sdm_warn <- function(msg, indent = 0L) {
  pad <- strrep("  ", indent)
  .sdm_msg(sprintf("%s  ! %s", pad, msg))
}

#' Start a named progress tracker and return a handle
#' @param name Label for this operation.
#' @param total Total number of steps (integer).
#' @return A list handle used with \code{sdm_progress_update} and \code{sdm_progress_done}.
#' @keywords internal
sdm_progress_start <- function(name, total = NULL) {
  list(name = name, n = 0L, total = total, start = proc.time()[["elapsed"]])
}

#' Advance a progress tracker by one step (prints nothing; step noted in handle)
#' @param handle Handle returned by \code{sdm_progress_start}.
#' @keywords internal
sdm_progress_update <- function(handle) {
  if (is.null(handle)) return(invisible(handle))
  handle$n <- handle$n + 1L
  invisible(handle)
}

#' Finish a progress tracker and print elapsed time
#' @param handle Handle returned by \code{sdm_progress_start}.
#' @keywords internal
sdm_progress_done <- function(handle) {
  if (is.null(handle) || !isTRUE(.alphasdm_env$verbose)) return(invisible(NULL))
  elapsed <- proc.time()[["elapsed"]] - handle$start
  .sdm_msg(sprintf("  \u2714 %s complete [%.1fs]", handle$name, elapsed))
  invisible(NULL)
}
