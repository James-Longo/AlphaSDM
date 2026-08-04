# format_data() reads the Alpha Earth coverage window from Earth Engine. These tests
# are about the data handling either side of that, so they seed the session cache the
# lookup writes to. That exercises the real code path with no network, and keeps the
# package free of a hardcoded window that would go stale.
local_coverage_window <- function(years = c(2017L, 2025L), env = parent.frame()) {
  old <- .alphasdm_env$year_range
  .alphasdm_env$year_range <- as.integer(years)
  withr::defer(.alphasdm_env$year_range <- old, envir = env)
  invisible(years)
}
