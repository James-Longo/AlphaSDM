# Package index

## Overview

- [`AlphaSDM`](https://james-longo.github.io/AlphaSDM/reference/AlphaSDM-package.md)
  [`AlphaSDM-package`](https://james-longo.github.io/AlphaSDM/reference/AlphaSDM-package.md)
  : AlphaSDM: Species Distribution Modeling using Alpha Earth Embeddings

## Connect to Earth Engine

- [`setup_gee()`](https://james-longo.github.io/AlphaSDM/reference/setup_gee.md)
  : Connect AlphaSDM to Google Earth Engine (one-time)
- [`gee_status()`](https://james-longo.github.io/AlphaSDM/reference/gee_status.md)
  : Report the Google Earth Engine connection status
- [`clear_gee_credentials()`](https://james-longo.github.io/AlphaSDM/reference/clear_gee_credentials.md)
  : Forget the saved Earth Engine project, and optionally sign out

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

## Earth Engine tasks and assets

- [`gee_tasks()`](https://james-longo.github.io/AlphaSDM/reference/gee_tasks.md)
  : Show recent AlphaSDM Earth Engine tasks
- [`gee_clean_assets()`](https://james-longo.github.io/AlphaSDM/reference/gee_clean_assets.md)
  : Remove leftover AlphaSDM temporary assets
