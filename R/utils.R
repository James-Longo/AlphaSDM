# Console output. Base R only, so the package pulls in no printing dependency.
# All output is message(), so suppressMessages() silences the package.

.alphasdm_env <- new.env(parent = emptyenv())

.sdm_msg <- function(text) message(text)

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
  if (is.null(handle)) return(invisible(NULL))
  elapsed <- proc.time()[["elapsed"]] - handle$start
  .sdm_msg(sprintf("  \u2714 %s complete [%.1fs]", handle$name, elapsed))
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
