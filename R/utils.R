# Console output and per-run state. Base R only, so the package pulls in no
# printing dependency.

# ---- Package-level state ----

.alphasdm_env <- new.env(parent = emptyenv())
.alphasdm_env$verbose <- TRUE
.alphasdm_env$standardization_active <- FALSE
.alphasdm_env$standardization_info_printed <- FALSE

#' Turn AlphaSDM console output on or off
#'
#' Every progress line the package prints is a message, not a warning or a print,
#' so this suppresses all of them at once. Errors are unaffected.
#'
#' @param verbose Logical. FALSE suppresses all non-error output.
#' @return `verbose`, invisibly. Called for its effect on package state.
#' @export
sdm_verbose <- function(verbose = TRUE) {
  .alphasdm_env$verbose <- isTRUE(verbose)
  invisible(verbose)
}

.sdm_msg <- function(text) {
  if (!isTRUE(.alphasdm_env$verbose)) return(invisible(NULL))
  message(text)
}

#' Print a section header
#' @param title Section title.
#' @return Nothing. Prints a message.
#' @noRd
sdm_section <- function(title) {
  rule <- strrep("\u2500", max(0, 60 - nchar(title) - 2))
  .sdm_msg(sprintf("\n\u250c\u2500 %s %s", title, rule))
}

#' Print a status line
#' @param msg Message text.
#' @param indent Indent level: 0 is top level, 1 is nested under a section.
#' @return Nothing. Prints a message.
#' @noRd
sdm_info <- function(msg, indent = 0L) {
  pad <- strrep("  ", indent)
  .sdm_msg(sprintf("%s  \u2023 %s", pad, msg))
}

#' Print a completion note
#' @param msg Message text.
#' @param indent Indent level: 0 is top level, 1 is nested under a section.
#' @return Nothing. Prints a message.
#' @noRd
sdm_done <- function(msg, indent = 0L) {
  pad <- strrep("  ", indent)
  .sdm_msg(sprintf("%s  \u2714 %s", pad, msg))
}

#' Print an advisory note
#'
#' Advisory only: this prints a message and does not signal a condition, so it
#' does not collect in `warnings()` and `suppressMessages()` hides it.
#'
#' @param msg Message text.
#' @param indent Indent level: 0 is top level, 1 is nested under a section.
#' @return Nothing. Prints a message.
#' @noRd
sdm_warn <- function(msg, indent = 0L) {
  pad <- strrep("  ", indent)
  .sdm_msg(sprintf("%s  ! %s", pad, msg))
}

#' Start a progress timer
#' @param name Label for the operation being timed.
#' @return A handle to pass to \code{sdm_progress_done}.
#' @noRd
sdm_progress_start <- function(name) {
  list(name = name, start = proc.time()[["elapsed"]])
}

#' Stop a progress timer and print how long it ran
#' @param handle Handle from \code{sdm_progress_start}.
#' @return Nothing. Prints a message.
#' @noRd
sdm_progress_done <- function(handle) {
  if (is.null(handle) || !isTRUE(.alphasdm_env$verbose)) return(invisible(NULL))
  elapsed <- proc.time()[["elapsed"]] - handle$start
  .sdm_msg(sprintf("  \u2714 %s complete [%.1fs]", handle$name, elapsed))
  invisible(NULL)
}

#' Reset the per-run console state
#'
#' `format_data()` prints its section header once per run, so the run tracks
#' whether it has printed yet. Register this with `on.exit()`. A run that aborts
#' would otherwise leave the flag set, and the next `format_data()` in that
#' session would print no header.
#' @noRd
reset_sdm_run_state <- function() {
  .alphasdm_env$standardization_active <- FALSE
  .alphasdm_env$standardization_info_printed <- FALSE
  invisible(NULL)
}

#' Close a run with a section header and the elapsed time
#' @param t_start Elapsed-time reading taken at the start of the run, in seconds.
#' @param title Section title to close with.
#' @return Nothing. Prints messages.
#' @noRd
sdm_finish <- function(t_start, title) {
  sdm_section(title)
  sdm_done(sprintf("Total elapsed time [%.1fs]", proc.time()[["elapsed"]] - t_start))
  .sdm_msg("")
}
