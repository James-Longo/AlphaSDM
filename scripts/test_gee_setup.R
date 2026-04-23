# AlphaSDM GEE Setup Integration Test
#
# This script tests the full clear → setup → verify cycle.
# It requires an interactive R session (browser auth will open).
#
# Usage:
#   source("scripts/test_gee_setup.R")

pkgload::load_all("/home/james-longo/Projects/AlphaSDM")

cat("\n=== Step 1: Clear all GEE credentials ===\n")
clear_gee_credentials()

# Verify cleared state
stopifnot(!file.exists(file.path(Sys.getenv("HOME"), ".config", "AlphaSDM", "config.json")))
stopifnot(is.null(getOption("AlphaSDM.gee_initialized")))
cat("  ✓ Credentials cleared successfully\n")

cat("\n=== Step 2: Run setup_gee() ===\n")
cat("  (This will open a browser for authentication)\n")
# Replace with your actual project ID
setup_gee(project = Sys.getenv("EARTHENGINE_PROJECT"))

# Verify setup state
stopifnot(file.exists(file.path(Sys.getenv("HOME"), ".config", "AlphaSDM", "config.json")))
stopifnot(isTRUE(getOption("AlphaSDM.gee_initialized")))
cat("  ✓ Setup completed successfully\n")

cat("\n=== Step 3: Verify GEE connection ===\n")
ee <- reticulate::import("ee")
info <- ee$Image(0)$getInfo()
stopifnot(!is.null(info))
cat("  ✓ GEE connection verified (ee$Image(0)$getInfo() returned successfully)\n")

cat("\n=== Step 4: Simulate new session (clear session cache) ===\n")
options(AlphaSDM.gee_initialized = NULL)
# This should auto-reconnect using saved project — no browser prompt
AlphaSDM:::ensure_gee_authenticated()
stopifnot(isTRUE(getOption("AlphaSDM.gee_initialized")))
cat("  ✓ Auto-reconnection works (simulated session restart)\n")

cat("\n=== All tests passed! ===\n")
