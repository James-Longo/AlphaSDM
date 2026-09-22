## Test environments

* Local: Ubuntu Linux, R 4.3.3
* win-builder: Windows, R-devel
* mac-builder: macOS (arm64), R-release

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* The words flagged as possibly misspelled in DESCRIPTION (AlphaEarth,
  embeddings, geospatial, et al.) are correct: AlphaEarth is the name of the
  Google DeepMind model whose embeddings the package uses.

## Notes for the reviewer

* AlphaSDM runs its computation on Google Earth Engine, which needs a free
  account and interactive sign-in. Examples that contact Earth Engine are
  therefore wrapped in `\dontrun{}`; the examples that run offline
  (`calculate_cbi()`, `calculate_classifier_metrics()`, `sdm_verbose()`) are
  not. Tests that need Earth Engine skip unless explicitly enabled.
* The vignette is precomputed from `vignettes/AlphaSDM.Rmd.orig`, so building
  it needs no credentials or network access.
* The Earth Engine Python client is declared with `reticulate::py_require()`;
  the package itself installs no software.
