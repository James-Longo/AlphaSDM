# AlphaSDM

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://james-longo.github.io/AlphaSDM/LICENSE.md)

AlphaSDM fits species distribution models and maps habitat suitability
at up to 10 m resolution, anywhere on Earth, from occurrence records
alone. It models species on the embeddings of AlphaEarth, Google
DeepMind’s geospatial foundation model, instead of environmental layers
you collect yourself, and runs every step on Google Earth Engine.

![Saguaro records and pseudo-absences on a Sentinel-2 image of
Tucson](reference/figures/README-points.png)![Saguaro habitat
suitability around Tucson at 30 m](reference/figures/README-map.png)

*Saguaro around Tucson, Arizona: GBIF records and pseudo-absences
(left), and the fitted habitat-suitability map at 30 m (right). The full
example is in
[`vignette("AlphaSDM")`](https://james-longo.github.io/AlphaSDM/articles/AlphaSDM.md).*

## Why AlphaSDM

- **No environmental layers.** There is nothing to find, download,
  reproject or align; the embeddings already describe every pixel.
- **Fine resolution everywhere.** 10 m pixels, every year from 2017, on
  land anywhere on Earth.
- **Nothing to download but the map.** Sampling, model fitting and
  prediction all run on Earth Engine, so a large study area costs your
  computer nothing.
- **An ensemble with calibration built in.** Support vector machine,
  random forest and boosted trees by default, scored with AUC, TSS and
  the Boyce index.
- **Explicit modelling choices.** Pseudo-absence placement follows
  Barbet-Massin et al. (2012), and AlphaSDM makes you choose the
  strategy rather than choosing it for you.

## AlphaEarth embeddings

[AlphaEarth Foundations](https://arxiv.org/abs/2507.22291) is a Google
DeepMind model that condenses optical, radar, lidar, climate and other
data into 64 numbers per 10 m pixel per year. The annual embeddings are
a public [Earth Engine
dataset](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_SATELLITE_EMBEDDING_V1_ANNUAL),
currently covering 2017 to 2025. Records are matched to the embeddings
for the year they were made.

## Installation

\
`# install.packages("pak")`\
`pak``::``pak``(``"James-Longo/AlphaSDM"``)`

## Earth Engine setup

AlphaSDM runs on your own Earth Engine account, which is free for
noncommercial use.

1.  [Register for Earth Engine](https://earthengine.google.com/signup/).
    This gives you a Cloud project ID.
2.  Connect once per machine. A browser window asks you to allow access,
    and the connection is remembered after that.

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`AlphaSDM`](https://github.com/James-Longo/AlphaSDM)`)`\
[`setup_gee`](https://james-longo.github.io/AlphaSDM/reference/setup_gee.md)`(``project ``=`` ``"your-project-id"``)`\
[`gee_status`](https://james-longo.github.io/AlphaSDM/reference/gee_status.md)`(``)``   ``# checks credentials, project and a live connection`

On a machine without a browser, use `setup_gee(auth_mode = "notebook")`
to paste a code instead.
[`clear_gee_credentials()`](https://james-longo.github.io/AlphaSDM/reference/clear_gee_credentials.md)
resets everything.

## Example

Download one year of saguaro records from GBIF, add pseudo-absences,
evaluate the default ensemble on a spatial holdout, and map suitability:

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`AlphaSDM`](https://github.com/James-Longo/AlphaSDM)`)`\
\
`url`` ``<-`` `[`paste0`](https://rdrr.io/r/base/paste.html)`(``"https://api.gbif.org/v1/occurrence/search?"``,`\
`              ``"scientificName=Carnegiea%20gigantea&year=2022"``,`\
`              ``"&hasCoordinate=true&hasGeospatialIssue=false"``,`\
`              ``"&coordinateUncertaintyInMeters=0,30"``,`\
`              ``"&decimalLongitude=-111.4,-110.6&decimalLatitude=31.9,32.6&limit=300"``)`\
`obs`` ``<-`` `[`do.call`](https://rdrr.io/r/base/do.call.html)`(``rbind``, `[`lapply`](https://rdrr.io/r/base/lapply.html)`(`[`c`](https://rdrr.io/r/base/c.html)`(``0``, ``300``)``, ``function``(``offset``)`\
`  ``jsonlite``::`[`fromJSON`](https://jeroen.r-universe.dev/jsonlite/reference/fromJSON.html)`(`[`paste0`](https://rdrr.io/r/base/paste.html)`(``url``, ``"&offset="``, ``offset``)``)``$``results``[`\
`    , `[`c`](https://rdrr.io/r/base/c.html)`(``"decimalLongitude"``, ``"decimalLatitude"``, ``"year"``)``]``)``)`\
\
`pres`` ``<-`` `[`format_data`](https://james-longo.github.io/AlphaSDM/reference/format_data.md)`(``obs``, coords ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"decimalLongitude"``, ``"decimalLatitude"``)``, year ``=`` ``"year"``)`\
`occ``  ``<-`` `[`generate_pseudo_absences`](https://james-longo.github.io/AlphaSDM/reference/generate_pseudo_absences.md)`(``pres``, aoi ``=`` ``"bbox"``, strategy ``=`` ``"combined"``,`\
`                                 n ``=`` `[`nrow`](https://rdrr.io/r/base/nrow.html)`(``pres``)``, aoi_year ``=`` ``2022``)`\
\
[`set.seed`](https://rdrr.io/r/base/Random.html)`(``1``)`\
`test`` ``<-`` ``stats``::`[`kmeans`](https://rdrr.io/r/stats/kmeans.html)`(``occ``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"longitude"``, ``"latitude"``)``]``, centers ``=`` ``5``)``$``cluster`` ``==`` ``1`\
`fit``  ``<-`` `[`evaluate_models`](https://james-longo.github.io/AlphaSDM/reference/evaluate_models.md)`(``occ``[``!``test``, ``]``, predict_coords ``=`` ``occ``[``test``, ``]``)`\
`fit``$``metrics``$``ensemble`\
\
`maps`` ``<-`` `[`generate_map`](https://james-longo.github.io/AlphaSDM/reference/generate_map.md)`(``occ``, aoi ``=`` ``"bbox"``, scale ``=`` ``30``, aoi_year ``=`` ``2022``,`\
`                     output_dir ``=`` ``"saguaro"``)`

[`generate_map()`](https://james-longo.github.io/AlphaSDM/reference/generate_map.md)
writes one GeoTIFF per model plus the ensemble. Maps download straight
from Earth Engine in tiles; a map Earth Engine will not compute that way
goes through its batch system and Google Drive instead, which is slower.

## Models

The default ensemble is `c("svm", "rf", "gbt")`. `methods =` also
accepts `"maxent"`, `"glm"`, `"cart"`, `"knn"`, `"mindist"` and
`"similarity"`, all fitted on Earth Engine; see
[`?evaluate_models`](https://james-longo.github.io/AlphaSDM/reference/evaluate_models.md).

## Related packages

biomod2, flexsdm, ENMeval, sdm and Wallace fit species distribution
models on environmental layers you supply; AlphaSDM replaces those
layers with one embedding dataset and moves the computation to Earth
Engine. blockCV builds spatial cross-validation folds, which pair well
with
[`evaluate_models()`](https://james-longo.github.io/AlphaSDM/reference/evaluate_models.md).
rgee is the general-purpose R interface to Earth Engine.

## Getting help

Report bugs and request features in [GitHub
issues](https://github.com/James-Longo/AlphaSDM/issues), or email
<james.longo.birds@gmail.com>. AlphaSDM is in active development, so
arguments and defaults may still change.

## Citation

Run `citation("AlphaSDM")` in R, or use GitHub’s “Cite this repository”
button, which reads
[`CITATION.cff`](https://james-longo.github.io/AlphaSDM/CITATION.cff).

## License

MIT; see
[LICENSE.md](https://james-longo.github.io/AlphaSDM/LICENSE.md). The
AlphaEarth embeddings are provided by Google under the terms of the
[Earth Engine
dataset](https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_SATELLITE_EMBEDDING_V1_ANNUAL).
