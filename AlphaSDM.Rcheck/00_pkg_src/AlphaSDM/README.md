# AlphaSDM: Making High-Resolution Species Distribution Modeling Faster, Easier, and More Accurate

`AlphaSDM` is an R package designed for large-scale, high-resolution species distribution modeling. It leverages **Google's Alpha Earth satellite embeddings** (64-dimensional dense vectors) to provide 10m-resolution mapping without the need for manual covariate selection or local data extraction.

All heavy computation—including embedding sampling, model training, and spatial prediction—is handled **server-side on Google Earth Engine (GEE)**, ensuring maximum performance even for millions of presence points.

**Google Earth Engine Access:**
`AlphaSDM` requires access to Google Earth Engine (GEE).
- **Free for Non-Commercial Use**: GEE is free for academic institutions, students, and nonprofit organizations for non-commercial research and educational purposes.
- **Registration**: You must register your Google account at [earthengine.google.com](https://earthengine.google.com/signup/).
- **Cloud Project**: You will need to create a **Google Cloud Project** and enable the **Earth Engine API**. When running `AlphaSDM` functions for the first time, use the `gee_project` argument to specify your Project ID.

---

## Key Features

*   **10m Resolution**: Native support for high-resolution 64-band Alpha Earth embeddings across the globe.
*   **Fully Server-Side**: No local covariate downloads. Data preparation and model execution happen on GEE's distributed infrastructure.
*   **Modeling Framework**: Support for diverse modeling approaches:
    *   **Classification**: Random Forest (RF), Gradient Boosted Trees (GBT), Support Vector Machines (SVM), Maxent.
    *   **Regression & Similarity**: Ridge Regression (Linear/Quadratic), Species Niche Centroid (Mean embedding dot-product).

---

## Installation

You can install `AlphaSDM` directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("James-Longo/AlphaSDM")
```

**Python Setup:**
`AlphaSDM` automatically handles Python dependencies via the `'reticulate'` package. It will check your active Python environment and install required packages (`earthengine-api`, `pandas`, `geopandas`, etc.) if they are missing.

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
