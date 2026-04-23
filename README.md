# AlphaSDM: Making High-Resolution Species Distribution Modeling Faster, Easier, and More Accurate

`AlphaSDM` is an R package designed for large-scale, high-resolution species distribution modeling. It leverages **Google's Alpha Earth satellite embeddings** (64-dimensional dense vectors) to provide 10m-resolution mapping without the need for manual covariate selection or local data extraction.

All heavy computation—including embedding sampling, model training, and spatial prediction—is handled **server-side on Google Earth Engine (GEE)**, ensuring maximum performance even for millions of presence points.

---

## Installation

You can install `AlphaSDM` directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("James-Longo/AlphaSDM")
```

---

## Google Earth Engine Setup

`AlphaSDM` requires a Google Earth Engine account. Setup is a one-time process:

### Prerequisites
- **GEE Account**: Register at [earthengine.google.com](https://earthengine.google.com/signup/). GEE is free for academic and non-commercial use.
- **Google Cloud Project**: Create a project at [console.cloud.google.com](https://console.cloud.google.com/) with the **Earth Engine API** enabled. Note your Project ID (e.g., `"my-ee-project"`).

### One-Time Setup

```r
library(AlphaSDM)

# First run only — installs Python dependencies and authenticates
setup_gee(project = "your-project-id")
```

On first run, `setup_gee()` will:
1. Install a Python environment with `earthengine-api` via `rgee` (you may need to restart R after this step)
2. Open a browser window for Google OAuth — just click **Allow**
3. Save your project ID locally so you never have to enter it again

**That's it.** All subsequent R sessions connect automatically — no further prompts or configuration needed.

### Resetting Credentials

To clear all saved credentials (e.g., for troubleshooting or switching accounts):

```r
clear_gee_credentials()

# Then re-run setup
setup_gee(project = "your-project-id")
```

---

## Key Features

*   **10m Resolution**: Native support for high-resolution 64-band Alpha Earth embeddings across the globe.
*   **Fully Server-Side**: No local covariate downloads. Data preparation and model execution happen on GEE's distributed infrastructure.
*   **Modeling Framework**: Support for diverse modeling approaches:
    *   **Classification**: Random Forest (RF), Gradient Boosted Trees (GBT), Support Vector Machines (SVM), Maxent.
    *   **Regression & Similarity**: Ridge Regression (Linear/Quadratic), Species Niche Centroid (Mean embedding dot-product).

---

## Quick Start

### 1. Format your data
Standardize your presence/absence records. Embeddings require a year for temporal alignment.

```r
library(AlphaSDM)

# Standardizes coordinates to (longitude, latitude, year, present)
# Coordinates must be provided as c(lon, lat)
formatted_data <- format_data(
  my_raw_df, 
  coords = c("lon", "lat"), 
  year = "obs_year", 
  presence = "occurrence"
)
```

### 2. Evaluate and Benchmark
Perform parallel evaluation at specific coordinates.

```r
metrics <- evaluate_models(
  data = train_data,
  predict_coords = test_data,
  scale = 10,
  aoi_year = 2024,
  methods = c("gbt", "maxent"),
  output_dir = "eval/run_1"
)

# Access metrics like AUC, TSS, and CBI
print(metrics$metrics$ensemble)
```

### 3. Generate Prediction Maps
Generate high-resolution suitability maps for a specific Area of Interest (AOI).

```r
# Define an AOI (Point + Radius or sf object)
aoi <- list(lat = 44.5, lon = -71.5, radius = 50000)

results <- generate_map(
  data = formatted_data,
  aoi = aoi,
  scale = 10,           
  aoi_year = 2023,      
  methods = c("rf", "ridge", "centroid"),
  output_dir = "results/my_species"
)
```

---

## License and Credits
This package leverages the **Alpha Earth Embedding** dataset provided by Google.
