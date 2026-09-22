# Package index

## Overview

- [`AlphaSDM`](https://james-longo.github.io/AlphaSDM/reference/AlphaSDM-package.md)
  [`AlphaSDM-package`](https://james-longo.github.io/AlphaSDM/reference/AlphaSDM-package.md)
  : AlphaSDM: Species Distribution Modeling using Alpha Earth Embeddings

## Connect to Earth Engine

- [`setup_gee()`](https://james-longo.github.io/AlphaSDM/reference/setup_gee.md)
  : Set Up Google Earth Engine for AlphaSDM (one-time)
- [`gee_status()`](https://james-longo.github.io/AlphaSDM/reference/gee_status.md)
  : Report the Google Earth Engine Connection Status
- [`clear_gee_credentials()`](https://james-longo.github.io/AlphaSDM/reference/clear_gee_credentials.md)
  : Clear All GEE Credentials and Configuration

## Prepare data

- [`format_data()`](https://james-longo.github.io/AlphaSDM/reference/format_data.md)
  : Standardize occurrence records for AlphaSDM
- [`generate_pseudo_absences()`](https://james-longo.github.io/AlphaSDM/reference/generate_pseudo_absences.md)
  : Generate pseudo-absences for presence-only data

## Model and map

- [`evaluate_models()`](https://james-longo.github.io/AlphaSDM/reference/evaluate_models.md)
  : Evaluate SDM models on Alpha Earth embeddings
- [`generate_map()`](https://james-longo.github.io/AlphaSDM/reference/generate_map.md)
  : Generate an SDM suitability map

## Evaluation metrics

- [`calculate_classifier_metrics()`](https://james-longo.github.io/AlphaSDM/reference/calculate_classifier_metrics.md)
  : Calculate discrimination and calibration metrics
- [`calculate_cbi()`](https://james-longo.github.io/AlphaSDM/reference/calculate_cbi.md)
  : Calculate the Continuous Boyce Index

## Earth Engine tasks and assets

- [`sdm_gee_status()`](https://james-longo.github.io/AlphaSDM/reference/sdm_gee_status.md)
  : Show the status of recent AlphaSDM Earth Engine batch tasks
- [`sdm_clean_assets()`](https://james-longo.github.io/AlphaSDM/reference/sdm_clean_assets.md)
  : Remove leftover AlphaSDM temporary assets
- [`sdm_verbose()`](https://james-longo.github.io/AlphaSDM/reference/sdm_verbose.md)
  : Turn AlphaSDM console output on or off
