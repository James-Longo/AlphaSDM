# AlphaSDM 0.1.1

## Breaking changes

* Presence-only data is no longer modelled with a silently generated
  background. `evaluate_models()` and `generate_map()` reject
  presence-only input; generate pseudo-absences explicitly first with
  the new `generate_pseudo_absences()`. `format_data()` still accepts
  presence-only records and tells you the next step.
* `evaluate_models()`/`generate_map()` lose the `count` argument; the
  background is whatever your data contains.

## New features

* `generate_pseudo_absences()`: strategies `random`, `disk`,
  `envelope`, `combined` following Barbet-Massin et al. (2012, Methods
  in Ecology and Evolution 3:327-338). You choose the strategy and AOI;
  the disk radius and environmental envelope threshold are estimated
  from your presences (and reported) unless you set them. Supply is
  guaranteed by redrawing, with the disk radius relaxing stepwise
  before any shortfall.
* `bg_replicates` (default TRUE): tree-family methods train on up to
  10 replicate balanced subsets of the absences and average their
  predictions, per the same paper's Table 1.
* Batch task waits now print the Earth Engine task-monitor links
  (Code Editor Tasks tab and the project's Cloud Console page) once per
  session.
* Failed prediction chunks are reported instead of silently returning
  NA.
