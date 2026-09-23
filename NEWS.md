# AlphaSDM (development version)

* Model settings are now given per model with `params`, using Earth Engine's
  argument names, for example `params = list(gbt = list(shrinkage = 0.01))`.
  This replaces the separate arguments `n_trees`, `shrinkage`, `svm_cost`,
  `knn_k` and the rest, which applied to several models at once and could
  silently override values equal to their defaults. Misspelt settings are now
  an error.
* Every model uses Earth Engine's defaults, except the tree counts Earth Engine
  requires (`rf` 500, `gbt` 150) and `knn`, which uses 15 neighbours because
  Earth Engine's single neighbour gives a two-value map. The SVM is now Earth
  Engine's default linear classifier, and boosted trees use Earth Engine's
  learning rate, which gives better-calibrated probabilities.
* Fixed classification SVMs (C_SVC, NU_SVC) returning inverted scores. Their
  probability depended on the order of the training rows, which Earth Engine
  does not preserve; the training table is now sorted before fitting.
* `generate_pseudo_absences()` now reads pseudo-absences from the same years as
  the presences, in the same proportions, instead of the latest embedding year.
  Pass `aoi_year` to place them all in one year.
* Models are trained on their rows in a fixed order, so retraining (which Earth
  Engine does for every map tile) gives the same model. Random forests are now
  exactly reproducible; boosted trees vary by at most about 0.005 inside Earth
  Engine.
* Removed naive Bayes, which discards negative values and so cannot use the
  signed embeddings.
* The vignette now tests models on the following year's records.

# AlphaSDM 0.2.0

First CRAN release.

* Maps download straight from Earth Engine in tiles sized to its 32 MB
  download limit, fetched two at a time and stitched into one GeoTIFF per
  model. Google Drive is used only for maps Earth Engine will not compute tile
  by tile.
* `generate_map()` computes the ensemble on Earth Engine and accepts
  `aoi = "bbox"`.
* The Earth Engine Python client is declared with `reticulate::py_require()`,
  so reticulate installs it when needed; AlphaSDM no longer installs Python
  itself and no longer depends on 'rgee'.
* The saved project ID lives in `tools::R_user_dir("AlphaSDM", "config")`.
  `clear_gee_credentials()` removes it and asks before signing out of Earth
  Engine.
* Connecting to Earth Engine no longer deletes old temporary assets; use
  `sdm_clean_assets()`.
* A getting-started vignette, `vignette("AlphaSDM")`, maps saguaro around
  Tucson from GBIF records.
* Removed unused code, including the remains of built-in cross-validation.
