## Resubmission

This is a resubmission. In this version I have:

* Single-quoted 'AlphaEarth' in the Title and Description, as requested.

Since the previous submission I have also fixed classification SVMs returning
inverted probability scores, made pseudo-absences use the same years as the
presences, and simplified the exported arguments and function names (see
NEWS.md). The package is still a new submission at version 0.2.0.

## Test environments

* Local: Ubuntu Linux, R 4.3.3
* win-builder: Windows, R-devel
* mac-builder: macOS (arm64), R-release

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.
* The other words flagged as possibly misspelled in DESCRIPTION (embeddings,
  geospatial, et al.) are spelled correctly.

## Notes for the reviewer

* AlphaSDM runs its computation on Google Earth Engine, which needs a free
  account and interactive sign-in. Examples that contact Earth Engine are
  therefore wrapped in `\dontrun{}`; the example that runs offline
  (`calculate_classifier_metrics()`) is not. Tests that need Earth Engine
  skip unless explicitly enabled.
* The vignette is precomputed from `vignettes/AlphaSDM.Rmd.orig`, so building
  it needs no credentials or network access.
* The Earth Engine Python client is declared with `reticulate::py_require()`;
  the package itself installs no software.
