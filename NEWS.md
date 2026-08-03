# AlphaSDM 0.1.0

First public release.

## Modeling

* `evaluate_models()` trains the ensemble on Google Earth Engine, runs spatial
  cross-validation and scores an independent set of coordinates, reporting
  AUC-ROC, AUC-PRG, TSS, balanced accuracy, the Continuous Boyce Index and
  Pearson correlation.
* `generate_map()` exports a continuous habitat-suitability GeoTIFF per model
  plus the ensemble, over a point-and-radius, a vector file or a pre-built
  `ee.Geometry`.
* Nine modeling methods are available through `methods =`, backed by
  `ee.Classifier` and `ee.Reducer`. The default ensemble is
  `c("svm", "rf", "gbt", "maxent")`.
* Spatial cross-validation defaults to `blockCV` spatial blocks, falling back to
  k-means coordinate clusters when `blockCV` is unavailable. `cv_method =
  "random"` selects non-spatial k-fold for interpolation-style tasks.
* Tree models train on a balanced 1:1 background by default (`balance_trees`),
  with `bg_ratio` available for an explicit absence:presence ratio.
* `weighted_ensemble = TRUE` weights ensemble members by cross-validated AUC.

## Earth Engine

* `setup_gee()` is a one-time, one-click connection: it reuses any Python
  environment that already provides `earthengine-api`, detects headless sessions
  and switches to the paste-a-code flow, and remembers the project ID.
* `gee_status()` reports the connection; `sdm_gee_status()` reports running
  server-side tasks; `clear_gee_credentials()` resets everything.
* Map export is Drive-free, tiled through `getDownloadURL` and mosaicked
  locally. Masked pixels (open water) round-trip through an explicit nodata
  sentinel so all-water tiles still return valid tiles. Set
  `ALPHASDM_TILE_CACHE` for resumable exports and `ALPHASDM_MAX_TILE_PX` to
  shrink tiles on restricted-quota projects.
* Point scoring samples embeddings at the evaluation coordinates and classifies
  the resulting FeatureCollection, which scales to high-abundance species.
  Chunks that hit a synchronous compute timeout retry through a batch export.
* Random forests and CART models can be persisted to a temporary Earth Engine
  asset so a whole-region export applies a stored model instead of retraining
  per tile. This falls back to the inline classifier if the batch scheduler is
  unavailable.
