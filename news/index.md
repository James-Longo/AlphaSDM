# Changelog

## AlphaSDM 0.2.0

First CRAN release.

- Maps download straight from Earth Engine in tiles sized to its 32 MB
  download limit, fetched two at a time and stitched into one GeoTIFF
  per model. Google Drive is used only for maps Earth Engine will not
  compute tile by tile.
- [`generate_map()`](https://james-longo.github.io/AlphaSDM/reference/generate_map.md)
  computes the ensemble on Earth Engine and accepts `aoi = "bbox"`.
- The Earth Engine Python client is declared with
  [`reticulate::py_require()`](https://rstudio.github.io/reticulate/reference/py_require.html),
  so reticulate installs it when needed; AlphaSDM no longer installs
  Python itself and no longer depends on ‘rgee’.
- The saved project ID lives in
  `tools::R_user_dir("AlphaSDM", "config")`.
  [`clear_gee_credentials()`](https://james-longo.github.io/AlphaSDM/reference/clear_gee_credentials.md)
  removes it and asks before signing out of Earth Engine.
- Connecting to Earth Engine no longer deletes old temporary assets; use
  [`sdm_clean_assets()`](https://james-longo.github.io/AlphaSDM/reference/sdm_clean_assets.md).
- A getting-started vignette,
  [`vignette("AlphaSDM")`](https://james-longo.github.io/AlphaSDM/articles/AlphaSDM.md),
  maps saguaro around Tucson from GBIF records.
- Removed unused code, including the remains of built-in
  cross-validation.
