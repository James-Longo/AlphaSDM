#' Split user settings into one list per method
#'
#' Shared by [evaluate_models()] and [generate_map()] so that the model you
#' evaluate is the model you map. Names that are not among `methods` are an
#' error, so a typo such as `gmb` does not silently leave a model on defaults.
#' @noRd
method_settings <- function(methods, params) {
  if (length(params) && (is.null(names(params)) || any(!nzchar(names(params)))))
    stop("`params` must be a named list with one entry per model, e.g. ",
         "list(gbt = list(shrinkage = 0.01)).", call. = FALSE)
  extra <- setdiff(names(params), methods)
  if (length(extra))
    stop(sprintf("`params` has settings for %s, which %s not in `methods`.",
                 paste(extra, collapse = ", "), if (length(extra) > 1) "are" else "is"),
         call. = FALSE)
  setNames(lapply(methods, function(m) params[[m]]), methods)
}

#' Reject presence-only training data with directions
#' @noRd
stop_if_presence_only <- function(data) {
  if (all(data$present == 1))
    stop("Presence-only data cannot be modelled: absence placement is a ",
         "modelling decision (Barbet-Massin et al. 2012). Generate ",
         "pseudo-absences first:\n  data <- generate_pseudo_absences(data, ",
         "aoi = ..., strategy = ...)\nSee ?generate_pseudo_absences for the ",
         "strategy recipes.", call. = FALSE)
}

#' Internal Unified GEE Training Pipeline
#' @noRd
fit_gee_models <- function(train_df, methods, scale, training_params,
                           bg_ratio = NULL, bg_replicates = FALSE, project = NULL) {
  ee <- reticulate::import("ee")

  # Optional class balancing: thin the background collection to a target
  # absence-to-presence ratio. Runs on Earth Engine, so no embeddings are downloaded.
  # This keeps the tree methods from collapsing onto the majority class, and keeps
  # the kNN vote fraction off the prevalence floor. It only ever removes points, so a
  # ratio looser than the data already is leaves the pool untouched.
  balance_bg <- function(bg_fc, n_pos, n_neg, seed = 42L) {
    if (is.null(bg_ratio) || n_pos <= 0L || n_neg <= 0L) return(bg_fc)
    target <- ceiling(as.numeric(bg_ratio) * n_pos)
    if (target >= n_neg) return(bg_fc)
    frac <- target / n_neg
    bg_fc$randomColumn("__bgsel", as.integer(seed))$filter(ee$Filter$lt("__bgsel", frac))
  }

  # Upload the points. Absence or pseudo-absence rows are required; both callers
  # reject presence-only input with stop_if_presence_only().
  sdm_section("Uploading training data to Google Earth Engine")
  pb_up <- sdm_progress_start("Uploading and sampling")

  # row_id travels with every training row so that training can sort by it:
  # Earth Engine does not preserve row order, and randomised models (boosted
  # trees, forests) otherwise differ slightly each time they are retrained,
  # which happens for every map tile. See train_gee_model().
  train_df$row_id <- seq_len(nrow(train_df))
  sdm_info(sprintf("Transferring %d coordinates ...", nrow(train_df)), indent = 1L)

  # Store the sampled training table once, in chunks of 5000 rows exported as
  # concurrent tasks. Everything downstream is lazy and re-evaluates its input
  # graph on every use, so training and every map tile then read stored values
  # instead of re-running the 64-band sampling.
  sdm_info("Computing the sampled training data (server-side) ...", indent = 1L)
  years  <- as.list(unique(as.integer(train_df$year)))
  chunks <- split(seq_len(nrow(train_df)), ceiling(seq_len(nrow(train_df)) / 5000))
  stored <- ee_store_tables(lapply(chunks, function(i) get_embeddings_at_fc(
    upload_points_to_gee(train_df[i, c("longitude", "latitude", "year", "present", "row_id")]),
    scale, properties = c("year", "present", "row_id"), geometries = TRUE, years = years)),
    project = project)
  sampled_fc     <- stored$fc
  training_asset <- stored$asset_ids
  pres_sampled   <- sampled_fc$filter(ee$Filter$eq("present", 1L))
  bg_all         <- sampled_fc$filter(ee$Filter$eq("present", 0L))
  n_pres         <- sum(train_df$present == 1)
  n_background   <- sum(train_df$present == 0)

  # Balanced-pool methods train on k replicate thinned subsets of the
  # absences and average (Barbet-Massin et al. 2012 Table 1: several runs
  # with few pseudo-absences; k from their measured asymptote). Replication
  # requires thinning to be active, since identical pools would be pointless.
  n_bal <- as.integer(ceiling(
    (if (is.null(bg_ratio)) 1 else as.numeric(bg_ratio)) * n_pres))
  k_rep <- if (isTRUE(bg_replicates) && !is.null(bg_ratio) &&
               n_bal < n_background)
    min(10L, as.integer(ceiling(10000 / max(n_bal, 1L)))) else 1L
  if (k_rep > 1L)
    sdm_info(sprintf(
      "Balanced pool: %d replicate subsets of %d absences (averaged)",
      k_rep, n_bal), indent = 1L)
  bal_pools <- lapply(seq_len(k_rep) - 1L, function(i)
    balance_bg(bg_all, n_pres, n_background, seed = 42L + i))
  bg_balanced <- bal_pools[[1]]

  sdm_progress_done(pb_up)

  # Earth Engine classifiers are lazy: clf$train() only builds a graph, so this loop
  # costs little. similarity is the exception and evaluates eagerly, in
  # train_gee_model().
  sdm_section(sprintf("Training %d model%s on Google Earth Engine",
                      length(methods), if (length(methods) == 1) "" else "s"))
  pb <- sdm_progress_start("Model training")
  models <- list()
  for (m in methods) {
    sdm_info(sprintf("Fitting %s ...", toupper(m)), indent = 1L)
    # Which pool a method trains on is a property of the method, declared in the
    # registry. See method_pool().
    fc_for_method <- switch(method_pool(m),
      presence = pres_sampled,
      balanced = pres_sampled$merge(bg_balanced),
      sampled_fc)
    rep_pools <- if (method_pool(m) == "balanced" && length(bal_pools) > 1L)
      bal_pools[-1] else list()
    # A large forest inline in the classify graph fails with "Computed value is
    # too large", so for large training sets the tree models Earth Engine can
    # store (rf, cart) are stored and loaded back.
    persist_m <- nrow(train_df) > 5000L && isTRUE(GEE_CLASSIFIER_METHODS[[m]]$persistable)
    models[[m]]   <- train_gee_model(fc_for_method, m, params = training_params[[m]],
                                     persist = persist_m, project = project)
    # No post-fit probe here, on purpose. Map export and scoring both classify on
    # their own and surface a malformed classifier anyway, and on a throttled tier
    # the extra synchronous round trip is the first thing to time out. That discarded
    # classifiers which had in fact trained.
    if (length(rep_pools) && isTRUE(models[[m]]$is_classifier) &&
        identical(models[[m]]$spec$transform, "none")) {
      models[[m]]$replicates <- lapply(rep_pools, function(bp)
        train_gee_model(pres_sampled$merge(bp), m,
                        params = training_params[[m]],
                        persist = FALSE, project = project)$trained)
    }
  }
  sdm_progress_done(pb)

  return(list(
    models   = models,
    metadata = list(
      n_presence       = n_pres,
      n_background     = n_background,
      methods          = methods,
      scale            = scale,
      # Temporary assets: stored classifiers plus the stored training table.
      # Removed by cleanup_classifier_assets() when the caller exits.
      classifier_assets = c(
        Filter(Negate(is.null), lapply(models, function(x) x$asset_id)),
        if (!is.null(training_asset)) as.list(training_asset)
      )
    )
  ))
}

#' Delete the temporary assets a fit created
#' @noRd
cleanup_classifier_assets <- function(train_res) {
  assets <- train_res$metadata$classifier_assets
  if (length(assets) > 0) for (a in assets) ee_delete_asset_quietly(a)
  invisible(NULL)
}


#' Generate an SDM suitability map
#'
#' Trains the model ensemble on Google Earth Engine and exports a continuous
#' habitat-suitability raster over an area of interest, one GeoTIFF per model plus
#' the ensemble. The maps download directly from Earth Engine in tiles; only a map
#' Earth Engine will not compute tile by tile goes through its batch system and
#' Google Drive, which is slower.
#'
#' @param data Data frame of training records with `longitude`, `latitude`, `year`
#'   and a `present` column (1 = presence; include 0 rows to supply real absences).
#' @param aoi Area of interest: a pre-built `ee.Geometry`, a list with `lon`/`lat`/`radius`,
#'   a path to a vector file readable by [sf::st_read()], or `"bbox"` for the
#'   bounding box of `data` (presences and absences).
#' @param scale Output resolution in metres (default 10).
#' @param output_dir Directory to write the GeoTIFF(s) to.
#' @param methods Character vector of models to ensemble. Defaults to
#'   `c("svm", "rf", "gbt")`; also accepts `maxent`, `glm` (logistic
#'   regression fitted server-side by IRLS with equal total class
#'   weights), `similarity`, `knn`, `cart`, `mindist`. MaxEnt and glm
#'   follow the regression-family recipe of Barbet-Massin et al. (2012):
#'   a large RANDOM pseudo-absence set suits them best (see
#'   `?generate_pseudo_absences`).
#' @param ensemble Logical; also export the ensemble mean map (default `TRUE`).
#' @param aoi_year Year of the embeddings to map. By default, the most common
#'   year in `data`.
#' @param bg_ratio Absence-to-presence ratio for the models that train on a
#'   balanced pool (rf, gbt and knn): their absences are thinned at random to
#'   `bg_ratio` times the number of presences. Default 1. `NULL` gives them
#'   every absence. svm and maxent always use every absence.
#' @param bg_replicates Logical (default TRUE). Train the balanced-pool
#'   methods (rf, gbt, knn) on k = min(10, ceil(10000/pool size))
#'   replicate thinned subsets of the absences and average their
#'   predictions (Barbet-Massin et al. 2012, Table 1: several runs when
#'   few pseudo-absences are used). Requires `bg_ratio` thinning to be
#'   active; methods on the full pool are never replicated.
#' @param params Named list of settings for individual models, using the
#'   argument names of the Earth Engine classifier, for example
#'   `list(gbt = list(shrinkage = 0.01), svm = list(kernelType = "RBF"))`.
#'   Every model uses Earth Engine's defaults except where Earth Engine needs a
#'   value or its default cannot work on these data: `rf` uses 500 trees and
#'   `gbt` 150, since Earth Engine requires a number, and `knn` uses 15
#'   neighbours, since Earth Engine's single neighbour gives a two-value map.
#'   See the Earth Engine reference for `ee.Classifier.smileRandomForest`,
#'   `smileGradientTreeBoost`, `libsvm`, `amnhMaxent`, `smileKNN` and
#'   `smileCart` for every available setting.
#' @param gee_project Optional Earth Engine project override (normally set via [setup_gee()]).
#' @return A named list of output file paths, with one `<method>_map` entry per
#'   model, plus `ensemble_map` when more than one method is requested.
#' @examples
#' \dontrun{
#' maps <- generate_map(occ, aoi = "bbox", scale = 30, output_dir = tempdir())
#' maps$ensemble_map
#' }
#' @export
generate_map <- function(data, aoi, scale = 10, output_dir = getwd(),
                         methods = c("svm", "rf", "gbt"), ensemble = TRUE,
                         aoi_year = NULL, bg_ratio = 1, bg_replicates = TRUE,
                         params = list(), gee_project = NULL) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  stop_if_presence_only(data)
  ensure_gee_authenticated(project = gee_project)
  t_total_start <- proc.time()[["elapsed"]]

  ee <- reticulate::import("ee")

  if (is.null(aoi_year)) aoi_year <- as.integer(names(which.max(table(data$year))))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  aoi_geom <- resolve_aoi(aoi, ee, data = data)

  method_params <- method_settings(methods, params)
  train_res <- fit_gee_models(data, methods, scale, method_params, bg_ratio = bg_ratio, bg_replicates = bg_replicates, project = gee_project)
  on.exit(cleanup_classifier_assets(train_res), add = TRUE)   # remove temp classifier assets on exit

  img_mosaic <- get_embedding_image(aoi_year)
  final_results <- list(methods = methods, model_metadata = train_res$metadata)
  want_ensemble <- isTRUE(ensemble) && length(methods) > 1L
  pb_map <- sdm_progress_start("Map generation")

  # One image with a band per model, plus the ensemble, so every map tile is a
  # single request that reads the embeddings once.
  #
  # The ensemble is the per-pixel mean of the members on their own scales, and
  # those scales do not agree. rf, gbt, maxent, knn and a classification svm
  # return probabilities on [0, 1]. similarity is a dot product against the
  # presence centroid, so it is signed and capped at that centroid's norm.
  # mindist is a difference of distances, signed and about twice as wide. A
  # regression svm (EPSILON_SVR) regresses the 0/1 label without clamping, so it
  # can fall outside [0, 1]. Mixing these lets the widest-spread member pull the mean around,
  # and the result is then not on a probability scale. Averaging within one
  # family, such as the default svm/rf/gbt tier, behaves.
  bands <- lapply(methods, function(m) predict_gee_map(train_res$models[[m]], img_mosaic)$rename(m))
  if (want_ensemble)
    bands <- c(bands, list(ee$Image$cat(bands)$reduce(ee$Reducer$mean())$rename("ensemble")))
  outputs <- c(methods, if (want_ensemble) "ensemble")
  paths   <- file.path(output_dir, paste0(outputs, ".tif"))
  sdm_info(sprintf("Exporting %s ...", paste(toupper(outputs), collapse = ", ")), indent = 1L)
  export_image(ee$Image$cat(bands), aoi_geom, scale, paths)
  for (i in seq_along(outputs)) final_results[[paste0(outputs[i], "_map")]] <- paths[i]

  sdm_progress_done(pb_map)
  
  sdm_finish(t_total_start, "Species processing finished")
  
  return(final_results)
}


#' Internal: score a coordinate data frame against trained models (FC-first, server-side)
#'
#' Samples the embeddings at the eval coordinates into a FeatureCollection, then classifies
#' that FC with each trained model. Classifying a finite point set this way is light and
#' scales to high-abundance species, unlike classifying the whole embedding image and
#' sampleRegions-ing it, whose per-pixel graph runs GEE out of memory for large models.
#' @noRd
predict_scores_internal <- function(predict_df, models, methods, scale,
                                    project = NULL, chunk_size = 4000L) {
  ee <- reticulate::import("ee")

  pred_cols <- paste0("pred_", methods)
  for (col in pred_cols) predict_df[[col]] <- NA_real_
  predict_df$.row_idx <- seq_len(nrow(predict_df))
  yrs <- as.list(sort(unique(as.integer(predict_df$year))))

  prop_num <- function(features, name) vapply(features, function(f) {
    v <- f$properties[[name]]
    if (is.null(v)) NA_real_ else as.numeric(v)
  }, numeric(1))
  fill_from_features <- function(features) {
    ridx <- prop_num(features, "row_idx")
    ok   <- !is.na(ridx)
    for (col in pred_cols) {
      vals <- prop_num(features, col)[ok]
      set  <- !is.na(vals)
      predict_df[[col]][ridx[ok][set]] <<- vals[set]
    }
  }

  # Sample the embeddings at the evaluation points first, then classify the resulting
  # collection. Classifying a few thousand feature vectors is light, and it scales to
  # species with many records. The alternative, classifying the whole embedding image
  # and sampling that, builds a per-pixel graph that runs Earth Engine out of memory
  # for a large model. Map export still takes the image path, in predict_gee_map();
  # this function only scores a finite set of points. A chunk that hits a genuine
  # compute timeout is retried through a batch export, which has a larger budget.
  score_features <- function(sub_df, batch) {
    chunk <- sub_df[, c("longitude", "latitude", "year", ".row_idx"), drop = FALSE]
    names(chunk)[names(chunk) == ".row_idx"] <- "row_idx"
    emb_fc <- get_embeddings_at_fc(upload_points_to_gee(chunk), scale,
                                   properties = c("year", "row_idx"),
                                   geometries = batch, years = yrs)  # exports need geometry
    scored <- predict_all_models_gee(emb_fc, models[methods])$select(as.list(c("row_idx", pred_cols)))
    if (batch) ee_read_via_batch(scored, project)$features else read_fc_paged(scored)$features
  }

  # One classify graph per chunk, read in pages so the 5000-feature cap does not
  # apply. A chunk that hits a compute timeout falls back to a batch export.
  chunk_size <- as.integer(chunk_size)
  n_rows     <- nrow(predict_df)
  for (i in seq(1L, n_rows, by = chunk_size)) {
    sub   <- predict_df[i:min(i + chunk_size - 1L, n_rows), , drop = FALSE]
    feats <- tryCatch(score_features(sub, batch = FALSE), error = function(e) e)
    if (inherits(feats, "error")) {
      if (is_gee_timeout(feats)) {
        sdm_warn("Scoring hit a GEE compute timeout; retrying that batch via export.", indent = 1L)
        feats <- tryCatch(score_features(sub, batch = TRUE), error = function(e) e)
        if (inherits(feats, "error")) { sdm_warn(sprintf("Batch scoring failed: %s", conditionMessage(feats)), indent = 1L); next }
      } else {
        # Never drop a chunk silently: an unlogged failure reads as NA
        # predictions with no explanation.
        sdm_warn(sprintf("Scoring chunk failed: %s",
                         substr(conditionMessage(feats), 1, 140)), indent = 1L)
        next
      }
    }
    fill_from_features(feats)
  }

  predict_df$.row_idx <- NULL
  return(predict_df)
}


#' Evaluate SDM models on Alpha Earth embeddings
#'
#' Trains the model ensemble on Google Earth Engine and scores an independent set
#' of coordinates. Training data must contain
#' absences: real ones, or pseudo-absences from [generate_pseudo_absences()].
#' Presence-only input is rejected with directions.
#'
#' @inheritParams generate_map
#' @param predict_coords Data frame of records to score, with `longitude`,
#'   `latitude` and `year` (as from [format_data()]). Include a `present` column
#'   to compute evaluation metrics.
#' @param scale Embedding resolution in metres (default 10, the native resolution).
#' @return A list containing `methods`, `model_metadata`,
#'   `point_predictions` and, when `predict_coords` has a `present`
#'   column, per-model and ensemble `metrics`. For cross-validation,
#'   split the data yourself and call this once per fold with the fold's
#'   holdout as `predict_coords`.
#' @examples
#' \dontrun{
#' # `occ` holds presences and absences, e.g. from generate_pseudo_absences().
#' test <- sample(nrow(occ), round(nrow(occ) / 5))
#' fit  <- evaluate_models(occ[-test, ], predict_coords = occ[test, ])
#' fit$metrics$ensemble
#' }
#' @export
evaluate_models <- function(data, predict_coords, scale = 10,
                            methods = c("svm", "rf", "gbt"), bg_ratio = 1,
                            bg_replicates = TRUE, params = list(), gee_project = NULL) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  stop_if_presence_only(data)
  need <- c("longitude", "latitude", "year")
  if (!all(need %in% names(predict_coords)))
    stop("`predict_coords` needs longitude, latitude and year columns; ",
         "format it with format_data().", call. = FALSE)
  ensure_gee_authenticated(project = gee_project)
  t_total_start <- proc.time()[["elapsed"]]

  method_params <- method_settings(methods, params)

  # Final model, trained on all the data.
  train_res <- fit_gee_models(data, methods, scale, method_params, bg_ratio = bg_ratio, bg_replicates = bg_replicates, project = gee_project)
  on.exit(cleanup_classifier_assets(train_res), add = TRUE)   # remove temp assets on any exit

  sdm_section(sprintf("Predicting at %d coordinates (server-side)", nrow(predict_coords)))
  pb_pred    <- sdm_progress_start("Prediction")
  final_pred <- predict_scores_internal(predict_coords, train_res$models, methods, scale, project = gee_project)
  sdm_progress_done(pb_pred)

  pred_cols <- paste0("pred_", methods)
  final_pred$pred_ensemble <- rowMeans(final_pred[, pred_cols, drop = FALSE], na.rm = TRUE)

  final_results <- list(
    methods           = c(methods, "ensemble"),
    metrics           = list(),
    point_predictions = final_pred,
    model_metadata    = train_res$metadata
  )

  if ("present" %in% names(final_pred)) {
    # Report evaluation points dropped for want of satellite coverage. Every model
    # reads the same image mask, so an NA in one model means NA in all of them.
    na_mask <- is.na(final_pred[[paste0("pred_", methods[1])]])
    if (any(na_mask)) {
      na_pres <- sum(na_mask & final_pred$present == 1)
      na_bg   <- sum(na_mask & final_pred$present == 0)
      sdm_warn(sprintf(
        "%d eval point%s dropped: no satellite coverage at prediction time (%d presence, %d background).",
        sum(na_mask), if (sum(na_mask) == 1) "" else "s", na_pres, na_bg
      ), indent = 1L)
    }

    for (m in c(methods, "ensemble")) {
      scores <- final_pred[[paste0("pred_", m)]]
      pos <- scores[final_pred$present == 1]
      neg <- scores[final_pred$present == 0]
      final_results$metrics[[m]] <- calculate_classifier_metrics(pos, neg)
    }
  }

  sdm_finish(t_total_start, "Species evaluation finished")

  return(final_results)
}
