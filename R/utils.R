#' Shared internal helper for live-updating timer
#' @keywords internal
log_time <- function(msg, start_time = NULL) {
    if (isTRUE(getOption("autoSDM.verbose_timer"))) {
        if (is.null(start_time)) {
             message("[autoSDM] ", msg, "...")
             return(Sys.time())
        } else {
             elapsed <- difftime(Sys.time(), start_time, units = "secs")
             message(sprintf("[autoSDM] %s... DONE (%0.2f seconds)", msg, as.numeric(elapsed)))
        }
    }
    return(NULL)
}
