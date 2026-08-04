# AlphaSDM

> [!IMPORTANT]
> AlphaSDM is in active development. Function arguments, defaults and outputs may
> still change between versions.
>
> If something doesn't work, or if it does, please email
> [james.longo@maine.edu](mailto:james.longo@maine.edu). Bugs and feature requests
> are also welcome as [GitHub issues](https://github.com/James-Longo/AlphaSDM/issues).

An R package for species distribution modeling at up to 10m resolution. AlphaSDM uses Google's [Alpha Earth satellite embeddings](https://arxiv.org/abs/2507.22291), 64-dimensional vectors that capture the environmental characteristics of any location on Earth, so you don't need to find, download, or align environmental layers yourself.

Everything runs on Google Earth Engine, from data extraction to model training to spatial prediction.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("James-Longo/AlphaSDM")
```

---

## Google Earth Engine Setup

AlphaSDM connects to Earth Engine with your personal Google account. You set
this up once per machine. After that, every R session connects automatically.

### Prerequisites (one-time, free)

1. Register for Earth Engine. Go to
   [earthengine.google.com/signup](https://earthengine.google.com/signup/) and
   sign in with your Google account. Earth Engine is free for noncommercial
   use in research, education, and nonprofit projects. Registration links your
   account to a Cloud project and gives you its Project ID (e.g.
   `"my-ee-project"`), which is the only thing you need to pass to AlphaSDM.

That's the only prerequisite. You do not need to create a service account,
manage a JSON key, or set up billing for standard noncommercial use.

### One-Time Setup

```r
library(AlphaSDM)

setup_gee(project = "your-project-id")
```

`setup_gee()` is designed to be run once and never thought about again:

1. Finds your Python automatically. If any Python on your machine already has
   `earthengine-api` (a conda env, a virtualenv, or a system Python with
   `pip install earthengine-api`), AlphaSDM uses it directly, with no downloads.
   It only builds a new environment as a last resort.
2. Authenticates with one browser click, no code to paste. On a desktop or
   laptop your browser opens, you click "Allow", and the credential is captured
   automatically. The saved credentials are long-lived; you are not asked again.
3. Remembers your project. The Project ID is saved locally.

Re-running `setup_gee()` is safe: if you are already connected it detects the
working credentials and returns immediately with *"Already connected to Earth
Engine. Nothing to do."*

### Check your connection

```r
gee_status()
#> ┌─ AlphaSDM: Google Earth Engine connection
#>     ‣ [OK  ] Python env : .../earthengine-api
#>     ‣ [OK  ] Credentials: present, user account (OAuth)
#>     ‣ [OK  ] Project    : my-ee-project
#>     ‣ [OK  ] Live check : connected
#>   ✔ Earth Engine is set up. No action needed.
```

### Headless / remote machines (SSH, HPC, containers)

On a machine with no local browser, `setup_gee()` detects this and prints a URL
to open on any device; you approve access and paste the short code back once. To
force this flow explicitly:

```r
setup_gee(project = "your-project-id", auth_mode = "notebook")
```

### Resetting Credentials

To switch accounts or troubleshoot:

```r
clear_gee_credentials()
setup_gee(project = "your-project-id")   # re-authenticate
```

---

## Key Features

*   Resolution: model habitat at up to 10m, anywhere on the globe, using Google's 64-band Alpha Earth satellite embeddings.
*   Fully server-side: no environmental data to download. All data extraction, model training, and prediction happens on Google Earth Engine.
*   Built-in models: nine modeling methods plus an ensemble, all trained and
    applied server-side. See [Built-in Models](#built-in-models) below.

---

## Built-in Models

Pass any of these to the `methods =` argument of `evaluate_models()` or
`generate_map()`. The default is `c("svm", "rf", "gbt", "maxent")`: four methods
that make different modelling assumptions, so averaging them is worth more than
averaging four variations on the same idea.

| `methods` value | Model | Earth Engine backend |
| --- | --- | --- |
| `"svm"` | Support Vector Machine (default EPSILON_SVR, RBF kernel) | [`ee.Classifier.libsvm`](https://developers.google.com/earth-engine/apidocs/ee-classifier-libsvm) |
| `"rf"` | Random Forest | [`ee.Classifier.smileRandomForest`](https://developers.google.com/earth-engine/apidocs/ee-classifier-smilerandomforest) |
| `"gbt"` | Gradient Boosted Trees | [`ee.Classifier.smileGradientTreeBoost`](https://developers.google.com/earth-engine/apidocs/ee-classifier-smilegradienttreeboost) |
| `"maxent"` | MaxEnt | [`ee.Classifier.amnhMaxent`](https://developers.google.com/earth-engine/apidocs/ee-classifier-amnhmaxent) |
| `"cart"` | Classification and Regression Tree | [`ee.Classifier.smileCart`](https://developers.google.com/earth-engine/apidocs/ee-classifier-smilecart) |
| `"knn"` | k-Nearest Neighbors | [`ee.Classifier.smileKNN`](https://developers.google.com/earth-engine/apidocs/ee-classifier-smileknn) |
| `"mindist"` | Minimum Distance to class centroid | [`ee.Classifier.minimumDistance`](https://developers.google.com/earth-engine/apidocs/ee-classifier-minimumdistance) |
| `"naivebayes"` | Naive Bayes (see note below) | [`ee.Classifier.smileNaiveBayes`](https://developers.google.com/earth-engine/apidocs/ee-classifier-smilenaivebayes) |
| `"similarity"` | Dot product against the mean presence embedding | `ee.Reducer.mean` |

Selecting more than one method also produces an ensemble map and score. It is
equal-weighted by default; `evaluate_models(weighted_ensemble = TRUE)` weights
members by their cross-validated AUC instead.

Two caveats worth knowing:

*   `"naivebayes"` is exposed for completeness but is not recommended.
    `smileNaiveBayes` assumes positive-integer features and discards negative
    inputs, so it collapses to roughly random performance on the signed Alpha
    Earth embeddings.
*   The lighter methods (`similarity`, `knn`, `cart`, `mindist`) are cheaper to fit
    but less expressive than the default four. `similarity` in particular is a single
    dot product against the mean presence embedding, with no fitted decision boundary.

---

## Quick Start

### 1. Format your data

Get your presence/absence records into the expected format. A year column is required so each record is matched to the embedding for the year it was recorded.

The Alpha Earth embeddings are annual and currently cover 2017 to 2025. `format_data()` drops any record dated outside that window and reports how many it removed. If the temporal gap does not matter for your question, set those records to 2017 to keep them and match them against the earliest available embedding.

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
  methods = c("gbt", "maxent")
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

## Citation

If AlphaSDM contributes to work you publish, please cite it. From R:

```r
citation("AlphaSDM")
```

The repository also carries a [`CITATION.cff`](CITATION.cff), so GitHub's
"Cite this repository" button produces BibTeX and APA entries directly.

## License and Credits

AlphaSDM is released under the [MIT License](LICENSE.md): you are free to use,
modify and redistribute it, including in commercial work, provided the copyright
notice is retained. If you use it in research, a citation is asked for as an
academic courtesy rather than a licence condition.

This package uses the [Alpha Earth Embedding](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_SATELLITE_EMBEDDING_V1_ANNUAL) dataset provided by Google.
