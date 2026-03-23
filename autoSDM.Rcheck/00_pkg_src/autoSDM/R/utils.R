#' Centralized logging with timestamps
#' @param msg Message string
#' @keywords internal
timestamp_message <- function(msg) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
}

#' Shared internal helper for live-updating timer
.run_command_with_timer <- function(label, active = FALSE) {
  if (active) timestamp_message(paste0("[autoSDM] ", label))
}
