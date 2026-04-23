# AlphaSDM

An R package for species distribution modeling at up to 10m resolution. AlphaSDM uses **Google's Alpha Earth satellite embeddings** — 64-dimensional vectors that capture the environmental characteristics of any location on Earth — so you don't need to find, download, or align environmental layers yourself.

Everything runs on **Google Earth Engine**, from data extraction to model training to spatial prediction.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("James-Longo/AlphaSDM")
```

---

## Google Earth Engine Setup

AlphaSDM needs a Google Earth Engine account. You only have to set this up once.

### Prerequisites
- **GEE Account**: Sign up at [earthengine.google.com](https://earthengine.google.com/signup/). Free for academic and non-commercial use.
- **Google Cloud Project**: Create one at [console.cloud.google.com](https://console.cloud.google.com/) and enable the **Earth Engine API**. Note your Project ID (e.g., `"my-ee-project"`).

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

**That's it.** Future R sessions connect automatically.

### Resetting Credentials

If you need to switch accounts or troubleshoot:

```r
clear_gee_credentials()

# Then re-run setup
setup_gee(project = "your-project-id")
```

---

## Key Features

*   **10m Resolution**: Model habitat at up to 10m resolution, anywhere on the globe, using Google's 64-band Alpha Earth satellite embeddings.
*   **Fully Server-Side**: No environmental data to download. All data extraction, model training, and prediction happens on Google Earth Engine.
*   **Built-in Models**: Five modeling methods ready to use:
    *   [**Random Forest (RF)**](https://developers.google.com/earth-engine/apidocs/ee-classifier-smilerandomforest)
    *   [**Gradient Boosted Trees (GBT)**](https://developers.google.com/earth-engine/apidocs/ee-classifier-smilegradienttreeboost)
    *   [**Maxent**](https://developers.google.com/earth-engine/apidocs/ee-classifier-amnhmaxent)
    *   [**Support Vector Machines (SVM)**](https://developers.google.com/earth-engine/apidocs/ee-classifier-libsvm)
    *   [**Similarity Search**](https://developers.google.com/earth-engine/apidocs/ee-reducer-mean) (dot product against the mean presence embedding)

---

## Quick Start

### 1. Format your data

Get your presence/absence records into the expected format. A year column is required for temporal alignment with the embeddings.

```r
library(AlphaSDM)

formatted_data <- format_data(
  my_raw_df, 
  coords = c("lon", "lat"), 
  year = "obs_year", 
  presence = "occurrence"
)
```

### 2. Evaluate models

Test model performance against a set of known coordinates:

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

### 3. Generate maps

Create maps for an area of interest. You can define the AOI in two ways:

```r
# Option 1: A center point with a radius (in meters)
aoi <- list(lat = 44.5, lon = -71.5, radius = 50000)

# Option 2: A path to any spatial file (Shapefile, GeoJSON, GeoPackage, KML, etc.)
aoi <- "path/to/my_study_area.shp"
```

```r
results <- generate_map(
  data = formatted_data,
  aoi = aoi,
  scale = 10,           
  aoi_year = 2023,      
  methods = c("rf", "similarity"),
  output_dir = "results/my_species"
)
```

---

## License and Credits

This package uses the [**Alpha Earth Embedding**](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_SATELLITE_EMBEDDING_V1_ANNUAL) dataset provided by Google.
