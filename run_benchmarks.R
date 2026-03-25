# AlphaSDM: Complete Benchmark Extraction (Final)
# 1. Setup Environment
env_name <- "r-AlphaSDM"
if (!reticulate::virtualenv_exists(env_name)) reticulate::virtualenv_create(env_name)
reticulate::use_virtualenv(env_name, required = TRUE)

if (!reticulate::py_module_available("ee")) {
  reticulate::py_install("earthengine-api", envname = env_name, pip = TRUE)
}

# 2. Source internal logic
if (file.exists("R/utils.R")) source("R/utils.R")
if (file.exists("R/gee_logic.R")) source("R/gee_logic.R")
source("scripts/extract_embeddings_to_drive.R")
library(data.table)

# 3. INITIALIZE EARTH ENGINE (Essential step)
message("AlphaSDM: Initializing Earth Engine...")
# Attempt to find Project ID from standard credentials
project_id <- NULL
creds_file <- "~/.config/earthengine/credentials"
if (file.exists(path.expand(creds_file))) {
  try({
    creds <- jsonlite::fromJSON(path.expand(creds_file))
    project_id <- creds$project
  }, silent = TRUE)
}
if (is.null(project_id)) project_id <- Sys.getenv("EARTHENGINE_PROJECT")

# Initialize
if (!is.null(project_id) && project_id != "") {
  rgee::ee_Initialize(project = project_id)
} else {
  rgee::ee_Initialize()
}

# 4. Load the 5.1M Row Dataset
print("AlphaSDM: Loading GeoPlant data (5.1M records)...")
df <- fread("../AlphaSDM_benchmarks/GeoPlant/PresenceOnlyOccurrences/PO_metadata_train.csv")
setnames(df, old=c("lat", "lon"), new=c("latitude", "longitude"))

# 5. Trigger the Year-by-Year extraction
print("AlphaSDM: Queuing safe (year-by-year) export tasks...")
tasks <- extract_embeddings_to_drive(df, "GeoPlant_Complete_Train", drive_folder="AlphaSDM_Benchmarking")

print("\nSuccess! Check GEE Task Manager.")
