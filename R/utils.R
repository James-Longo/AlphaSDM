#' Shared internal helper for live-updating timer
#' @importFrom processx process
.run_command_with_timer <- function(command, args, label, active = FALSE) {
  if (!active || !requireNamespace("processx", quietly = TRUE)) {
    return(system2(command, args = args))
  }
  
  start_time <- Sys.time()
  p <- processx::process$new(command, args, stdout = "|", stderr = "|")
  
  while (p$is_alive()) {
    elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")))
    cat("\r[autoSDM] ", label, " running for ", elapsed, " seconds...", sep = "")
    Sys.sleep(1)
  }
  
  # Check if it failed
  status <- p$get_exit_status()
  
  cat("\r[autoSDM] ", label, " completed in ", 
      round(as.numeric(difftime(Sys.time(), start_time, units = "secs"))), 
      " seconds.                 \n", sep = "")
  
  return(status)
}
