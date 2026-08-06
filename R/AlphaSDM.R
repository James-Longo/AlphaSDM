#' Build the per-method hyperparameter list
#'
#' Shared by [evaluate_models()] and [generate_map()] so that the model you
#' cross-validate is the model you map. Each method starts from the shared tree
#' arguments and then takes its own overrides.
#'
#' Every override is conditional on the corresponding argument still holding its
#' default, so an explicit value from the caller always wins.
#' @noRd
build_method_params <- function(methods, n_trees, min_leaf_population, bag_fraction,
                                shrinkage, max_nodes, variables_per_split,
                                svm_type, svm_kernel, svm_cost, svm_gamma,
                                maxent_beta, maxent_features,
                                knn_k, knn_search_method, knn_metric) {
  base <- list(numberOfTrees = n_trees, minLeafPopulation = min_leaf_population,
               bagFraction = bag_fraction, shrinkage = shrinkage,
               maxNodes = max_nodes, variablesPerSplit = variables_per_split)
  p <- setNames(lapply(methods, function(m) base), methods)
  has <- function(m) m %in% methods

  # gbt: many shallow trees taking small steps.
  if (has("gbt") && n_trees == 100L)    p$gbt$numberOfTrees <- 150L
  if (has("gbt") && shrinkage == 0.005) p$gbt$shrinkage     <- 0.05

  # rf: the opposite. It averages many low-bias trees, so it wants depth the shared
  # defaults do not give, and a split subset near sqrt of the 64 embedding bands.
  if (has("rf") && n_trees == 100L)                 p$rf$numberOfTrees      <- 500L
  if (has("rf") && is.null(variables_per_split))    p$rf$variablesPerSplit  <- 8L
  if (has("rf") && max_nodes == 6L)                 p$rf$maxNodes           <- NULL
  if (has("rf") && min_leaf_population == 5L)       p$rf$minLeafPopulation  <- 1L

  if (has("svm")) {
    p$svm$svmType    <- svm_type
    p$svm$kernelType <- svm_kernel
    p$svm$cost       <- svm_cost
    p$svm$gamma      <- svm_gamma
  }
  # MaxEnt, in the style of ENMeval: a regularisation multiplier and feature classes.
  if (has("maxent")) p$maxent <- modifyList(p$maxent, maxent_tuning_params(maxent_beta, maxent_features))

  # k sets the resolution of the kNN surface, not just its smoothness. In PROBABILITY
  # mode the score is the positive-vote fraction, so it can take only k + 1 distinct
  # values. Raise it alongside bg_ratio: k should stay a small share of the balanced
  # pool, or the vote smooths away the signal.
  if (has("knn") && !is.null(knn_k)) p$knn$k <- as.integer(knn_k)
  # searchMethod and metric are unset by default, so Earth Engine chooses them. Both
  # matter. The documentation warns that results vary by search method "for distance
  # ties and probability values", and KD_TREE ignores the metric, as does AUTO at low
  # dimensions. Set searchMethod to LINEAR_SEARCH or COVER_TREE whenever the metric
  # is meant to apply.
  if (has("knn") && !is.null(knn_search_method)) p$knn$searchMethod <- as.character(knn_search_method)
  if (has("knn") && !is.null(knn_metric))        p$knn$metric       <- as.character(knn_metric)

  p
}

#' The default model ensemble
#'
#' Four methods that make different modelling assumptions, so their errors are only
#' partly correlated and averaging them is worth more than averaging four variations
#' on one idea. The lighter methods (similarity, knn, cart, mindist) stay available
#' via `methods=`.
#' @noRd
DEFAULT_METHODS <- c("svm", "rf", "gbt", "maxent")

#' Internal Unified GEE Training Pipeline
#' @noRd
fit_gee_models <- function(train_df, methods, aoi_geom, scale, aoi_year, training_params, count = NULL,
                           bg_ratio = NULL, persist_classifier = TRUE, project = NULL) {
  # persist_classifier defaults to TRUE. Which methods actually store a model is a
  # property of the method, held in the registry's `persistable` field, so a caller
  # never has to reason about it: svm, gbt and maxent ignore the argument. Storing a
  # whole-region map's tree model once lets every export tile apply the stored forest
  # instead of retraining it, which is the largest export saving on a throttled tier.
  # Cross-validation passes FALSE, since a model stored per fold is pure overhead.
  ee <- reticulate::import("ee")

  # Optional class balancing: thin the background collection to a target
  # absence-to-presence ratio. Runs on Earth Engine, so no embeddings are downloaded.
  # This keeps the tree methods from collapsing onto the majority class, and keeps
  # the kNN vote fraction off the prevalence floor. It only ever removes points, so a
  # ratio looser than the data already is leaves the pool untouched.
  balance_bg <- function(bg_fc, n_pos, n_neg) {
    if (is.null(bg_ratio) || n_pos <= 0L || n_neg <= 0L) return(bg_fc)
    target <- ceiling(as.numeric(bg_ratio) * n_pos)
    if (target >= n_neg) return(bg_fc)
    frac <- target / n_neg
    bg_fc$randomColumn("__bgsel", 42L)$filter(ee$Filter$lt("__bgsel", frac))
  }

  # Upload the points and, for presence-only input, generate the background. The
  # background is created and kept on Earth Engine.
  sdm_section("Uploading training data to Google Earth Engine")
  pb_up <- sdm_progress_start("Uploading and sampling")

  if (all(train_df$present == 1)) {
    pres_df  <- train_df[, c("longitude", "latitude", "year", "present")]
    n_pres   <- nrow(pres_df)
    bg_count <- if (is.null(count)) n_pres else as.integer(count)

    lons <- pres_df$longitude; lats <- pres_df$latitude
    sdm_info(sprintf("Coordinate extent: Lon [%.2f, %.2f]  Lat [%.2f, %.2f]",
                     min(lons), max(lons), min(lats), max(lats)), indent = 1L)
    sdm_info(sprintf("Transferring %d presence coordinates ...", n_pres), indent = 1L)

    all_years   <- as.list(unique(c(as.integer(pres_df$year), as.integer(aoi_year))))
    presence_fc <- upload_points_to_gee(pres_df)

    pres_sampled <- get_embeddings_at_fc(presence_fc, scale,
                                         properties = c("year", "present"),
                                         years      = all_years)

    # Draw the background with randomPoints on Earth Engine. The graph stays the same
    # size whatever n is, which avoids the 10 MB payload limit that uploading sampled
    # embeddings would hit.
    n_bg <- min(n_pres, bg_count)
    sdm_info(sprintf("Sampling %d background points (server-side) ...", n_bg), indent = 1L)
    bg_fc <- generate_background_fc_gee(aoi_year, n_bg, aoi_geom)

    bg_props   <- c("year", "present")
    bg_sampled <- get_embeddings_at_fc(bg_fc, scale,
                                       properties = bg_props, years = list(as.integer(aoi_year)))

    bg_balanced <- balance_bg(bg_sampled, n_pres, n_bg)
    sampled_fc  <- pres_sampled$merge(bg_sampled)
    n_background <- n_bg
  } else {
    upload_df    <- train_df[, c("longitude", "latitude", "year", "present")]
    sdm_info(sprintf("Transferring %d coordinates ...", nrow(upload_df)), indent = 1L)

    upload_fc    <- upload_points_to_gee(upload_df)
    sampled_fc   <- get_embeddings_at_fc(upload_fc, scale,
                                         properties = c("year", "present"),
                                         years      = as.list(unique(as.integer(upload_df$year))))
    pres_sampled <- sampled_fc$filter(ee$Filter$eq("present", 1L))
    # Supplied absences serve as the background pool, as in the branch above.
    n_pres       <- sum(train_df$present == 1)
    n_background <- sum(train_df$present == 0)
    bg_balanced  <- balance_bg(sampled_fc$filter(ee$Filter$eq("present", 0L)), n_pres, n_background)
  }

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
    models[[m]]   <- train_gee_model(fc_for_method, m, params = training_params[[m]],
                                     persist = persist_classifier, project = project)
    # No post-fit probe here, on purpose. Map export and scoring both classify on
    # their own and surface a malformed classifier anyway, and on a throttled tier
    # the extra synchronous round trip is the first thing to time out. That discarded
    # classifiers which had in fact trained.
  }
  sdm_progress_done(pb)

  return(list(
    models   = models,
    metadata = list(
      n_presence       = n_pres,
      n_background     = n_background,
      methods          = methods,
      scale            = scale,
      # Temporary classifier assets written by persist_classifier, empty if none.
      # Remove them with cleanup_classifier_assets() once prediction is finished.
      classifier_assets = Filter(Negate(is.null), lapply(models, function(x) x$asset_id))
    )
  ))
}

#' Delete temporary classifier assets created by persist_classifier
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
#' the ensemble. Export is Drive-free (tiled `getDownloadURL`).
#'
#' @param data Data frame of training records with `longitude`, `latitude`, `year`
#'   and a `present` column (1 = presence; include 0 rows to supply real absences).
#' @param aoi Area of interest: a pre-built `ee.Geometry`, a list with `lon`/`lat`/`radius`,
#'   or a path to a vector file readable by [sf::st_read()].
#' @param scale Output resolution in metres (default 10).
#' @param output_dir Directory to write the GeoTIFF(s) to.
#' @param methods Character vector of models to ensemble. Defaults to
#'   `c("svm", "rf", "gbt", "maxent")`; also accepts `similarity`, `knn`, `cart`, `mindist`.
#' @param ensemble Logical; also export the ensemble mean map (default `TRUE`).
#' @param aoi_year Year of the Alpha Earth mosaic to sample (default 2023).
#' @param count Number of background points to generate for presence-only data.
#' @param bg_ratio Optional absence:presence ratio for the balanced background pool
#'   (the methods whose registry entry declares `pool = "balanced"`). Overrides
#'   `balance_trees` when set.
#' @param balance_trees Logical (default `TRUE`). When `TRUE`, rf/gbt and knn train on a
#'   balanced 1:1 background while svm/maxent use the full background; `FALSE` gives
#'   the trees all background points.
#' @param n_trees,min_leaf_population,bag_fraction,shrinkage,max_nodes,variables_per_split
#'   Tree-model (rf/gbt) hyperparameters.
#' @param svm_type,svm_kernel,svm_cost,svm_gamma libsvm hyperparameters (default
#'   EPSILON_SVR / RBF / cost 10 / gamma 0.05).
#' @param maxent_beta,maxent_features MaxEnt regularisation multiplier and feature classes.
#' @param knn_k Neighbours for kNN (default 15). Also fixes the output resolution:
#'   the surface can take only `k + 1` distinct values. Raise alongside `bg_ratio`.
#' @param knn_search_method kNN neighbour search: `"AUTO"`, `"LINEAR_SEARCH"`,
#'   `"KD_TREE"` or `"COVER_TREE"`. Note `KD_TREE` ignores `knn_metric`.
#' @param knn_metric kNN distance metric: `"EUCLIDEAN"`, `"MAHALANOBIS"`,
#'   `"MANHATTAN"` or `"BRAYCURTIS"`. Only honoured for search methods that use it.
#' @param persist_classifier Logical; whether to persist internally-persistable classifiers
#'   (currently RF/CART) to a temporary GEE asset so a whole-region export applies a stored
#'   model instead of retraining the forest on every tile. Defaults to `TRUE`; the package
#'   ignores it for methods that cannot persist (SVM/GBT/MaxEnt), so most callers should
#'   leave it unset.
#' @param gee_project Optional Earth Engine project override (normally set via [setup_gee()]).
#' @return A named list of output file paths, with one `<method>_map` entry per
#'   model, plus `ensemble_map` when more than one method is requested.
#' @export
generate_map <- function(data, aoi, scale = 10, output_dir = getwd(),
                         methods = NULL, ensemble = TRUE, aoi_year = NULL, count = NULL, bg_ratio = NULL,
                         balance_trees = TRUE,
                         n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                         shrinkage = 0.005, max_nodes = 6L, variables_per_split = NULL,
                         svm_type = "EPSILON_SVR", svm_kernel = "RBF", svm_cost = 10, svm_gamma = 0.05,
                         maxent_beta = 1, maxent_features = "auto",
                         knn_k = NULL, knn_search_method = NULL, knn_metric = NULL,
                         persist_classifier = TRUE,
                         gee_project = NULL) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)
  on.exit(reset_sdm_run_state(), add = TRUE)
  t_total_start <- proc.time()[["elapsed"]]

  ee <- reticulate::import("ee")

  if (is.null(aoi_year)) aoi_year <- 2023
  if (is.null(methods)) methods <- DEFAULT_METHODS

  # Balanced 1:1 background by default. Set balance_trees = FALSE for the full
  # background, or pass a numeric bg_ratio to override both. See evaluate_models().
  if (is.null(bg_ratio) && isTRUE(balance_trees)) bg_ratio <- 1
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  aoi_geom <- NULL
  if (inherits(aoi, "python.builtin.object")) {
    aoi_geom <- aoi                                   # pre-built ee.Geometry
  } else if (is.list(aoi) && !is.null(aoi$lat)) {
    aoi_geom <- ee$Geometry$Point(c(as.numeric(aoi$lon), as.numeric(aoi$lat)))$buffer(as.numeric(aoi$radius))
  } else if (is.character(aoi) && file.exists(aoi)) {
    aoi_sf <- sf::st_read(aoi, quiet = TRUE)
    aoi_geom <- rgee::sf_as_ee(aoi_sf)$geometry()
  } else {
    stop("`aoi` must be an ee.Geometry, a list with lon/lat/radius, or a path to a ",
         "readable vector file. Got: ", paste(class(aoi), collapse = "/"), call. = FALSE)
  }

  method_params <- build_method_params(
    methods, n_trees, min_leaf_population, bag_fraction, shrinkage, max_nodes,
    variables_per_split, svm_type, svm_kernel, svm_cost, svm_gamma,
    maxent_beta, maxent_features, knn_k, knn_search_method, knn_metric)
  train_res <- fit_gee_models(data, methods, aoi_geom, scale, aoi_year, method_params, count = count, bg_ratio = bg_ratio, persist_classifier = persist_classifier, project = gee_project)
  on.exit(cleanup_classifier_assets(train_res), add = TRUE)   # remove temp classifier assets on exit

  img_mosaic <- get_embedding_image(aoi_year)
  final_results <- list(methods = methods, model_metadata = train_res$metadata)
  want_ensemble <- isTRUE(ensemble) && length(methods) > 1L
  pb_map <- sdm_progress_start("Map generation")
  member_tifs <- character(0)
  for (m in methods) {
    sdm_info(sprintf("Exporting %s ...", toupper(m)), indent = 1L)
    pred_img <- predict_gee_map(train_res$models[[m]], img_mosaic)
    tif_path <- file.path(output_dir, paste0(m, ".tif"))
    # Export without Google Drive. The export gives back a single raster when the
    # region is one tile and a directory of tiles otherwise, so record what it
    # actually produced rather than the name it was asked for.
    produced <- tryCatch(
      export_image_tiled(pred_img, aoi_geom, scale, tif_path),
      error = function(e) {
        sdm_warn(sprintf("Export of %s failed: %s", m, conditionMessage(e))); NULL
      }
    )
    if (!is.null(produced)) {
      final_results[[paste0(m, "_map")]] <- produced
      member_tifs <- c(member_tifs, produced)
    }
  }

  # Ensemble map: the per-pixel mean of the members. The member rasters were just
  # written to output_dir on the same grid, so average them here rather than ask
  # Earth Engine to classify every pixel with every model a second time.
  #
  # The members are averaged on their own scales, and those scales do not agree.
  # rf, gbt, maxent and knn return probabilities on [0, 1]. similarity is a dot
  # product against the presence centroid, so it is signed and capped at that
  # centroid's norm. mindist is a difference of distances, signed and about twice as
  # wide. svm under its EPSILON_SVR default regresses the 0/1 label without clamping,
  # so it can fall outside [0, 1] altogether. Mixing these lets the widest-spread
  # member pull the mean around, and the result is then not on a probability scale.
  # Averaging within one family, such as the default svm/rf/gbt/maxent tier, behaves.
  if (want_ensemble) {
    written <- Filter(function(p) file.exists(p), member_tifs)
    if (length(written) < 2L) {
      sdm_warn("Ensemble skipped: fewer than two member maps were exported.", indent = 1L)
    } else if (all(!dir.exists(written))) {
      # Single-raster members: average them directly.
      sdm_info("Writing ENSEMBLE ...", indent = 1L)
      tif_path <- file.path(output_dir, "ensemble.tif")
      ok <- tryCatch({
        layers <- lapply(written, function(p) stars::read_stars(p, proxy = FALSE))
        acc <- layers[[1]]
        stacked <- simplify2array(lapply(layers, function(r) r[[1]]))
        acc[[1]][] <- apply(stacked, seq_len(length(dim(stacked)) - 1L), mean, na.rm = TRUE)
        stars::write_stars(acc, tif_path); TRUE
      }, error = function(e) { sdm_warn(sprintf("Ensemble failed: %s", conditionMessage(e))); FALSE })
      if (ok) final_results$ensemble_map <- tif_path
    } else {
      # Tiled members: average tile by tile, so the whole region is never in memory
      # at once. Only tiles present for every member are combined.
      sdm_info("Writing ENSEMBLE tiles ...", indent = 1L)
      ens_dir <- file.path(output_dir, "ensemble_tiles")
      dir.create(ens_dir, showWarnings = FALSE, recursive = TRUE)
      common <- Reduce(intersect, lapply(written, function(d) basename(list.files(d, "\\.tif$"))))
      n_ok <- 0L
      for (nm in common) {
        ok <- tryCatch({
          layers <- lapply(written, function(d) stars::read_stars(file.path(d, nm), proxy = FALSE))
          acc <- layers[[1]]
          stacked <- simplify2array(lapply(layers, function(r) r[[1]]))
          acc[[1]][] <- apply(stacked, seq_len(length(dim(stacked)) - 1L), mean, na.rm = TRUE)
          stars::write_stars(acc, file.path(ens_dir, nm)); TRUE
        }, error = function(e) FALSE)
        if (isTRUE(ok)) n_ok <- n_ok + 1L
      }
      if (n_ok > 0L) {
        final_results$ensemble_map <- ens_dir
        sdm_info(sprintf("%d ensemble tiles written.", n_ok), indent = 2L)
      } else {
        sdm_warn("Ensemble tiles could not be written.", indent = 1L)
      }
    }
  }

  sdm_progress_done(pb_map)
  
  sdm_finish(t_total_start, "Species processing finished")
  
  return(final_results)
}

#' Score a pre-sampled embedding FeatureCollection
#'
#' Cross-validation scores the same validation background in every fold, and only
#' the classifier changes between folds. Sampling those embeddings once and then
#' classifying the stored collection per fold turns `cv_folds` full 64-band
#' sampleRegions passes into one. Rows are matched back by the 0-based `row_idx`
#' carried through Earth Engine.
#' @noRd
score_sampled_fc <- function(emb_fc, models, methods, n_rows, async = FALSE, project = NULL) {
  pred_cols <- paste0("pred_", methods)
  out <- as.data.frame(matrix(NA_real_, nrow = n_rows, ncol = length(methods)))
  names(out) <- pred_cols

  scored <- predict_all_models_gee(emb_fc, models[methods])$select(as.list(c("row_idx", pred_cols)))
  feats  <- if (isTRUE(async)) ee_table_to_info_async(scored, project)$features else read_fc_paged(scored)$features

  for (f in feats) {
    ridx <- suppressWarnings(as.integer(f$properties[["row_idx"]]))
    if (length(ridx) != 1L || is.na(ridx) || ridx < 0L || ridx >= n_rows) next
    for (m in methods) {
      v <- f$properties[[paste0("pred_", m)]]
      if (!is.null(v)) out[ridx + 1L, paste0("pred_", m)] <- as.numeric(v)
    }
  }
  out
}

#' Internal: score a coordinate data frame against trained models (FC-first, server-side)
#'
#' Samples the embeddings at the eval coordinates into a FeatureCollection, then classifies
#' that FC with each trained model. Classifying a finite point set this way is light and
#' scales to high-abundance species, unlike classifying the whole embedding image and
#' sampleRegions-ing it, whose per-pixel graph runs GEE out of memory for large models.
#' @noRd
predict_scores_internal <- function(predict_df, models, methods, scale, aoi_year,
                                    async = FALSE, project = NULL, chunk_size = 4000L) {
  ee <- reticulate::import("ee")
  if (!"year" %in% names(predict_df)) predict_df$year <- aoi_year

  pred_cols <- paste0("pred_", methods)
  for (col in pred_cols) predict_df[[col]] <- NA_real_
  predict_df$.row_idx <- seq_len(nrow(predict_df))
  yrs <- as.list(sort(unique(as.integer(predict_df$year))))

  fill_from_features <- function(features) {
    for (f in features) {
      ridx <- as.integer(f$properties[["row_idx"]])
      if (is.null(ridx) || is.na(ridx)) next
      for (m in methods) {
        val <- f$properties[[paste0("pred_", m)]]
        if (!is.null(val)) predict_df[ridx, paste0("pred_", m)] <<- as.numeric(val)
      }
    }
  }

  # Sample the embeddings at the evaluation points first, then classify the resulting
  # collection. Classifying a few thousand feature vectors is light, and it scales to
  # species with many records. The alternative, classifying the whole embedding image
  # and sampling that, builds a per-pixel graph that runs Earth Engine out of memory
  # for a large model. Map export still takes the image path, in predict_gee_map();
  # this function only scores a finite set of points. A chunk that hits a genuine
  # compute timeout is retried through a batch export, which has a larger budget.
  score_features <- function(sub_df, use_async) {
    chunk <- sub_df[, c("longitude", "latitude", "year", ".row_idx"), drop = FALSE]
    names(chunk)[names(chunk) == ".row_idx"] <- "row_idx"
    emb_fc <- get_embeddings_at_fc(upload_points_to_gee(chunk), scale,
                                   properties = c("year", "row_idx"),
                                   geometries = use_async, years = yrs)  # geometry needed for export
    scored <- predict_all_models_gee(emb_fc, models[methods])$select(as.list(c("row_idx", pred_cols)))
    if (use_async) ee_table_to_info_async(scored, project)$features else read_fc_paged(scored)$features
  }

  # One classify graph per chunk, read in pages so the 5000-feature cap does not
  # apply. A chunk that hits a compute timeout falls back to a batch export.
  chunk_size <- as.integer(chunk_size)
  n_rows     <- nrow(predict_df)
  for (i in seq(1L, n_rows, by = chunk_size)) {
    sub   <- predict_df[i:min(i + chunk_size - 1L, n_rows), , drop = FALSE]
    feats <- tryCatch(score_features(sub, use_async = isTRUE(async)), error = function(e) e)
    if (inherits(feats, "error")) {
      if (is_gee_timeout(feats)) {
        sdm_warn("Scoring hit a GEE compute timeout; retrying that batch via export.", indent = 1L)
        feats <- tryCatch(score_features(sub, use_async = TRUE), error = function(e) e)
        if (inherits(feats, "error")) { sdm_warn(sprintf("Batch scoring failed: %s", conditionMessage(feats)), indent = 1L); next }
      } else next
    }
    fill_from_features(feats)
  }

  predict_df$.row_idx <- NULL
  return(predict_df)
}


#' Evaluate SDM models on Alpha Earth embeddings
#'
#' Trains the model ensemble on Google Earth Engine and either cross-validates it
#' and/or scores an independent set of coordinates. Training data may be
#' presence-only (background is generated automatically) or presence/absence
#' (the supplied absences are used directly).
#'
#' @param data Data frame of training records with `longitude`, `latitude`, `year`
#'   and a `present` column (1 = presence; include 0 rows to supply real absences).
#' @param predict_coords Optional data frame of coordinates to score. Include a
#'   `present` column to compute evaluation metrics on it.
#' @param scale Embedding resolution in metres (default 10, the native resolution).
#' @param methods Character vector of models to ensemble. Defaults to
#'   `c("svm", "rf", "gbt", "maxent")`; also accepts `similarity`, `knn`,
#'   `cart`, `mindist`.
#' @param aoi_year Year of the Alpha Earth mosaic to sample (default 2023).
#' @param count Number of background points to generate for presence-only data.
#' @param bg_ratio Optional absence:presence ratio; downsamples the balanced background
#'   pool (the methods whose registry entry declares `pool = "balanced"`) to
#'   `bg_ratio * n_presence`. Overrides `balance_trees` when set.
#' @param balance_trees Logical (default `TRUE`). When `TRUE`, rf/gbt and knn train on a
#'   balanced 1:1 background while svm/maxent use the full background; set `FALSE`
#'   to train the trees on all background points too.
#' @param n_trees,min_leaf_population,bag_fraction,shrinkage,max_nodes,variables_per_split
#'   Tree-model (rf/gbt) hyperparameters.
#' @param svm_type,svm_kernel,svm_cost,svm_gamma libsvm hyperparameters (default
#'   EPSILON_SVR / RBF / cost 10 / gamma 0.05).
#' @param maxent_beta,maxent_features MaxEnt regularisation multiplier and feature
#'   classes (`"auto"` or a combination of L/Q/H/P/T).
#' @param knn_k Neighbours for kNN (default 15). Also fixes the output resolution:
#'   the surface can take only `k + 1` distinct values. Raise alongside `bg_ratio`.
#' @param knn_search_method kNN neighbour search: `"AUTO"`, `"LINEAR_SEARCH"`,
#'   `"KD_TREE"` or `"COVER_TREE"`. Note `KD_TREE` ignores `knn_metric`.
#' @param knn_metric kNN distance metric: `"EUCLIDEAN"`, `"MAHALANOBIS"`,
#'   `"MANHATTAN"` or `"BRAYCURTIS"`. Only honoured for search methods that use it.
#' @param cv_folds Number of spatial cross-validation folds (0 to skip CV).
#' @param cv_method Fold assignment: `"block"` (spatial blocks), `"kmeans"`
#'   (coordinate clusters), or `"random"` (non-spatial k-fold).
#' @param block_size Optional block size (metres) for `cv_method = "block"`.
#' @param weighted_ensemble Logical; if `TRUE`, weight the ensemble by CV AUC
#'   (requires `cv_folds > 0`). Default `FALSE` (equal weights).
#' @param async Logical; use asynchronous GEE export for large prediction sets.
#' @param persist_classifier Logical; persist internally-persistable classifiers (RF/CART)
#'   to a temporary GEE asset. Defaults to `FALSE` here: cross-validation refits per fold,
#'   where per-fold persistence is pure overhead. (Map generation defaults it `TRUE`.)
#' @param gee_project Optional Earth Engine project override (normally set via [setup_gee()]).
#' @param options Named list of advanced options; `batch_size` sets how many
#'   coordinates are scored per Earth Engine request (default 4000).
#' @return A list containing `methods`, `model_metadata`, `cv_metrics` and
#'   `ensemble_weights`; when `predict_coords` is supplied it also includes
#'   `point_predictions` and (if a `present` column is present) per-model and
#'   ensemble `metrics`.
#' @export
evaluate_models <- function(data, predict_coords = NULL, scale = 10,
                            methods = NULL, aoi_year = NULL, count = NULL, bg_ratio = NULL,
                            balance_trees = TRUE,
                            n_trees = 100L, min_leaf_population = 5L, bag_fraction = 0.5,
                            shrinkage = 0.005, max_nodes = 6L, variables_per_split = NULL,
                            svm_type = "EPSILON_SVR", svm_kernel = "RBF", svm_cost = 10, svm_gamma = 0.05,
                            maxent_beta = 1, maxent_features = "auto",
                            knn_k = NULL, knn_search_method = NULL, knn_metric = NULL,
                            cv_folds = 5L, cv_method = "block", block_size = NULL,
                            weighted_ensemble = FALSE, async = FALSE,
                            persist_classifier = FALSE,
                            gee_project = NULL,
                            options = list()) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)
  on.exit(reset_sdm_run_state(), add = TRUE)
  t_total_start <- proc.time()[["elapsed"]]

  ee <- reticulate::import("ee")

  if (is.null(aoi_year)) aoi_year <- 2023
  if (is.null(methods)) methods <- DEFAULT_METHODS

  # Tree models train on a balanced 1:1 background by default, which matters most on
  # imbalanced presence and background data. svm and maxent always see the full
  # background. Set balance_trees = FALSE to give the trees all of it, or pass a
  # numeric bg_ratio, an absence-to-presence ratio, to override both.
  if (is.null(bg_ratio) && isTRUE(balance_trees)) bg_ratio <- 1

  # Points scored per Earth Engine request. Lower it when a large job trips the
  # synchronous compute budget partway through a chunk.
  score_chunk <- if (!is.null(options$batch_size)) as.integer(options$batch_size) else 4000L

  cv_folds <- as.integer(cv_folds)
  if (weighted_ensemble && cv_folds == 0L) {
    cv_folds <- 5L
    sdm_warn("weighted_ensemble = TRUE requires cross-validation; cv_folds overridden to 5.")
  }
  if (is.null(predict_coords) && cv_folds == 0L) {
    stop("Either predict_coords or cv_folds > 0 must be provided.")
  }

  method_params <- build_method_params(
    methods, n_trees, min_leaf_population, bag_fraction, shrinkage, max_nodes,
    variables_per_split, svm_type, svm_kernel, svm_cost, svm_gamma,
    maxent_beta, maxent_features, knn_k, knn_search_method, knn_metric)
  # Area of interest: the bounding box of the prediction targets, or of the training
  # data when no prediction coordinates were given.
  ref_df   <- if (!is.null(predict_coords)) predict_coords else data
  bbox     <- c(min(ref_df$longitude), min(ref_df$latitude),
                max(ref_df$longitude), max(ref_df$latitude))
  aoi_geom <- ee$Geometry$Rectangle(bbox)

  cv_metrics       <- NULL
  ensemble_weights <- NULL

  if (cv_folds > 0L) {
    sdm_section(sprintf("Cross-validation (%d spatial folds)", cv_folds))

    pres_df <- data[data$present == 1, , drop = FALSE]
    n_pres  <- nrow(pres_df)

    if (cv_folds > n_pres) {
      stop(sprintf(
        "cv_folds (%d) exceeds the number of presence points (%d). Reduce cv_folds or provide more training data.",
        cv_folds, n_pres
      ))
    }

    # Spatial folds from blockCV by default. Its blocks control spatial
    # autocorrelation better than the contiguous regions k-means produces. Falls back
    # to k-means when blockCV is not installed.
    pres_df$cv_fold <- assign_cv_folds(pres_df, cv_folds, method = cv_method, block_size = block_size)
    sdm_info(sprintf("CV folds: %s (%d folds)",
                     switch(cv_method, block = "blockCV blocks", kmeans = "k-means clusters",
                            random = "random k-fold", cv_method), cv_folds), indent = 1L)

    # Generate the validation background and download its coordinates once, before
    # the fold loop.
    val_bg_n  <- if (is.null(count)) n_pres else as.integer(count)
    bg_fc_gee <- generate_background_fc_gee(aoi_year, val_bg_n, aoi_geom)
    bg_info   <- retry_curl_download(bg_fc_gee$getInfo())
    val_bg_df <- do.call(rbind, lapply(bg_info$features, function(f) {
      data.frame(longitude = f$geometry$coordinates[[1]],
                 latitude  = f$geometry$coordinates[[2]],
                 year      = aoi_year)
    }))

    # The validation background is the same in every fold, and only the classifier
    # changes, so sample its embeddings once and reuse the stored collection. Falls
    # back to sampling per fold when the batch scheduler is unavailable, which it can
    # be on a throttled tier.
    val_bg_df$row_idx <- seq_len(nrow(val_bg_df)) - 1L
    bg_emb <- get_embeddings_at_fc(
      upload_points_to_gee(val_bg_df[, c("longitude", "latitude", "year", "row_idx")]),
      scale, properties = c("year", "row_idx"), geometries = TRUE,
      years = list(as.integer(aoi_year)))
    bg_mat <- tryCatch(ee_materialize_fc_async(bg_emb, project = gee_project, max_minutes = 20),
                       error = function(e) {
                         sdm_warn(sprintf("Could not cache the validation background (%s); sampling it per fold.",
                                          conditionMessage(e)), indent = 1L)
                         NULL
                       })
    if (!is.null(bg_mat)) on.exit(ee_delete_asset_quietly(bg_mat$asset_id), add = TRUE)

    cv_fold_aucs <- matrix(NA_real_, nrow = cv_folds, ncol = length(methods),
                           dimnames = list(NULL, methods))

    for (k in seq_len(cv_folds)) {
      sdm_info(sprintf("Fold %d / %d", k, cv_folds), indent = 1L)

      train_k <- pres_df[pres_df$cv_fold != k, c("longitude", "latitude", "year", "present"), drop = FALSE]
      test_k  <- pres_df[pres_df$cv_fold == k, c("longitude", "latitude", "year"),             drop = FALSE]

      res_k        <- fit_gee_models(train_k, methods, aoi_geom, scale, aoi_year, method_params, count = count, bg_ratio = bg_ratio, persist_classifier = persist_classifier, project = gee_project)
      pres_scored  <- predict_scores_internal(test_k,    res_k$models, methods, scale, aoi_year, async = async, project = gee_project, chunk_size = score_chunk)
      bg_scored    <- if (!is.null(bg_mat)) {
        score_sampled_fc(bg_mat$fc, res_k$models, methods, nrow(val_bg_df), async = async, project = gee_project)
      } else {
        predict_scores_internal(val_bg_df, res_k$models, methods, scale, aoi_year, async = async, project = gee_project, chunk_size = score_chunk)
      }
      cleanup_classifier_assets(res_k)   # per-fold temp classifier assets

      na_pres_cv <- is.na(pres_scored[[paste0("pred_", methods[1])]])
      na_bg_cv   <- is.na(bg_scored[[paste0("pred_",  methods[1])]])
      if (any(na_pres_cv) || any(na_bg_cv)) {
        sdm_warn(sprintf(
          "CV fold %d: %d presence and %d background points dropped: no satellite coverage.",
          k, sum(na_pres_cv), sum(na_bg_cv)
        ), indent = 2L)
      }
      for (m in methods) {
        pos <- pres_scored[[paste0("pred_", m)]]; pos <- pos[!is.na(pos)]
        neg <- bg_scored[[paste0("pred_", m)]];   neg <- neg[!is.na(neg)]
        if (length(pos) > 0 && length(neg) > 0)
          cv_fold_aucs[k, m] <- calculate_classifier_metrics(pos, neg)$auc_roc
      }
    }

    cv_aucs    <- colMeans(cv_fold_aucs, na.rm = TRUE)
    cv_metrics <- as.list(cv_aucs)
    sdm_info(paste("CV AUC:", paste(sprintf("%s:%.3f", names(cv_aucs), cv_aucs), collapse = "  ")))

    if (weighted_ensemble) {
      raw_wts <- pmax(cv_aucs - 0.5, 0)
      if (sum(raw_wts) == 0) {
        sdm_warn("All CV AUCs <= 0.5; falling back to equal weights")
        ensemble_weights <- setNames(rep(1 / length(methods), length(methods)), methods)
      } else {
        ensemble_weights <- raw_wts / sum(raw_wts)
      }
      sdm_info(paste("Ensemble weights:", paste(sprintf("%s:%.3f", names(ensemble_weights), ensemble_weights), collapse = "  ")))
    }
  }

  # Final model, trained on all the data.
  train_res <- fit_gee_models(data, methods, aoi_geom, scale, aoi_year, method_params, count = count, bg_ratio = bg_ratio, persist_classifier = persist_classifier, project = gee_project)
  on.exit(cleanup_classifier_assets(train_res), add = TRUE)   # remove temp assets on any exit

  # Cross-validation only: return without scoring any prediction coordinates.
  if (is.null(predict_coords)) {
    sdm_finish(t_total_start, "Species evaluation finished")
    return(list(
      methods          = c(methods, "ensemble"),
      model_metadata   = train_res$metadata,
      cv_metrics       = cv_metrics,
      ensemble_weights = ensemble_weights
    ))
  }

  sdm_section(sprintf("Predicting at %d coordinates (server-side)", nrow(predict_coords)))
  pb_pred    <- sdm_progress_start("Prediction")
  final_pred <- predict_scores_internal(predict_coords, train_res$models, methods, scale, aoi_year, async = async, project = gee_project, chunk_size = score_chunk)
  sdm_progress_done(pb_pred)

  pred_cols <- paste0("pred_", methods)
  if (weighted_ensemble && !is.null(ensemble_weights)) {
    wts <- ensemble_weights[methods]
    final_pred$pred_ensemble <- as.numeric(as.matrix(final_pred[, pred_cols, drop = FALSE]) %*% wts)
  } else {
    final_pred$pred_ensemble <- rowMeans(final_pred[, pred_cols, drop = FALSE], na.rm = TRUE)
  }

  final_results <- list(
    methods           = c(methods, "ensemble"),
    metrics           = list(),
    point_predictions = final_pred,
    model_metadata    = train_res$metadata,
    cv_metrics        = cv_metrics,
    ensemble_weights  = ensemble_weights
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
      dropped_coords <- final_pred[na_mask, c("longitude", "latitude", "present"), drop = FALSE]
      dropped_coords$class <- ifelse(dropped_coords$present == 1, "presence", "background")
      for (j in seq_len(nrow(dropped_coords))) {
        sdm_info(sprintf("  Dropped %s: lon=%.4f, lat=%.4f",
                         dropped_coords$class[j],
                         dropped_coords$longitude[j],
                         dropped_coords$latitude[j]), indent = 2L)
      }
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
