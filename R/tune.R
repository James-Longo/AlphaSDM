# =============================================================================
#  Hyperparameter / method sweep with spatial-CV model selection.
#  tune_models() trains a grid of method x hyperparameter configurations under
#  spatial cross-validation (blockCV by default), pools out-of-fold predictions,
#  and ranks them by AUC/TSS — the package analogue of an ENMeval/caret grid
#  search. Embeddings are sampled ONCE (server-side) and reused across all
#  configurations, so the whole sweep costs one embedding extraction.
# =============================================================================

#' Cartesian product of candidate parameter values -> list of parameter lists
#'
#' Each argument is a vector or list of candidate values (use a list to include
#' NULL, e.g. `maxNodes = list(NULL, 50L)` for unlimited vs capped depth).
#' @keywords internal
.cfg_grid <- function(...) {
  args <- list(...)
  if (length(args) == 0) return(list(list()))
  idx <- expand.grid(lapply(args, seq_along), KEEP.OUT.ATTRS = FALSE)
  lapply(seq_len(nrow(idx)), function(i)
    setNames(lapply(names(args), function(nm) args[[nm]][[idx[i, nm]]]), names(args)))
}

#' Compact human-readable label for a configuration
#' @keywords internal
.cfg_label <- function(cfg) {
  if (length(cfg) == 0) return("default")
  paste(vapply(names(cfg), function(nm) {
    v <- cfg[[nm]]
    sprintf("%s=%s", nm, if (is.null(v)) "inf" else as.character(v))
  }, character(1)), collapse = ",")
}

#' Default per-method tuning grids
#'
#' Sensible, non-exhaustive grids covering the levers that matter for Alpha Earth
#' embeddings: SVM type/kernel/cost/gamma, RF depth/trees/mtry, GBT trees/shrinkage/
#' depth, and MaxEnt regularization x feature classes. Override via `grids`.
#' @keywords internal
default_tuning_grids <- function(methods) {
  g <- list()
  if ("svm" %in% methods) g$svm <- c(
    .cfg_grid(svmType = "EPSILON_SVR", kernelType = c("RBF", "SIGMOID"), cost = c(10, 100), gamma = c(0.05, 0.5)),
    .cfg_grid(svmType = "EPSILON_SVR", kernelType = "LINEAR", cost = c(10, 100)),
    .cfg_grid(svmType = "NU_SVR",      kernelType = "RBF", cost = 10, gamma = c(0.05, 0.5)),
    .cfg_grid(svmType = "C_SVC",       kernelType = "RBF", cost = 100, gamma = 0.5))
  if ("rf" %in% methods) g$rf <- .cfg_grid(
    numberOfTrees = c(250L, 500L), maxNodes = list(NULL, 50L),
    variablesPerSplit = c(8L, 16L), minLeafPopulation = 1L)
  if ("gbt" %in% methods) g$gbt <- .cfg_grid(
    numberOfTrees = c(100L, 150L), shrinkage = c(0.005, 0.05), maxNodes = c(6L, 20L))
  if ("maxent" %in% methods) g$maxent <- c(
    list(maxent_tuning_params(1, "auto")),
    do.call(c, lapply(c("LQ", "LQH", "LQHP"),
                      function(fc) lapply(c(0.5, 1, 2, 4), function(b) maxent_tuning_params(b, fc)))))
  if ("knn"  %in% methods) g$knn  <- .cfg_grid(k = c(5L, 15L, 30L, 50L))
  if ("cart" %in% methods) g$cart <- .cfg_grid(maxNodes = list(NULL, 50L), minLeafPopulation = c(1L, 5L))
  if ("mindist" %in% methods) g$mindist <- .cfg_grid(metric = c("euclidean", "cosine"))
  if ("similarity" %in% methods) g$similarity <- list(list())   # niche centroid (no hyperparameters)
  # Output-mode checks: the registry's default mode is continuous (good for AUC), but
  # let the sweep confirm it where an alternative is plausible (rf/gbt PROBABILITY vs
  # REGRESSION, C_SVC PROBABILITY+invert vs MULTIPROBABILITY).
  if ("rf"  %in% methods) g$rf  <- c(g$rf,  .cfg_grid(numberOfTrees = 500L, minLeafPopulation = 1L, output = "REGRESSION"))
  if ("gbt" %in% methods) g$gbt <- c(g$gbt, .cfg_grid(numberOfTrees = 150L, shrinkage = 0.05, output = "REGRESSION"))
  if ("svm" %in% methods) g$svm <- c(g$svm, .cfg_grid(svmType = "C_SVC", kernelType = "RBF", cost = 100, gamma = 0.5, output = "MULTIPROBABILITY"))
  # naivebayes is intentionally omitted: GEE smileNaiveBayes discards negative inputs,
  # so it collapses to ~0.5 on the [-1,1] Alpha Earth embeddings.
  g
}

#' Score a test FeatureCollection with a trained model, returned as named scores
#' @keywords internal
score_test_fc <- function(model_res, te_fc, id = "rid") {
  ee   <- reticulate::import("ee")
  spec <- if (!is.null(model_res$spec)) model_res$spec else GEE_CLASSIFIER_METHODS[[model_res$method]]
  info <- retry_curl_download(
    te_fc$classify(model_res$trained)$select(list(id, spec$score))$getInfo())
  ids <- vapply(info$features, function(f) as.character(as.integer(f$properties[[id]])), character(1))
  read1 <- if (identical(spec$transform, "mindist_raw")) {
    # minimumDistance RAW: score is array [d_absence, d_presence] -> suitability = d_abs - d_pres
    function(f) { v <- f$properties[[spec$score]]; as.numeric(v[[1]]) - as.numeric(v[[2]]) }
  } else if (identical(spec$transform, "multiprob_last")) {
    # MULTIPROBABILITY: score is an array of per-class probs; presence is the last element
    function(f) { v <- f$properties[[spec$score]]; as.numeric(v[[length(v)]]) }
  } else {
    function(f) as.numeric(f$properties[[spec$score]])
  }
  sc <- setNames(vapply(info$features, read1, numeric(1)), ids)
  if (identical(spec$transform, "invert")) sc <- 1 - sc
  sc
}

#' Score a test FeatureCollection with a trained reducer (similarity = centroid dot-product)
#' @keywords internal
score_similarity_fc <- function(model_res, te_fc, id = "rid") {
  ee <- reticulate::import("ee")
  emb_cols <- sprintf("A%02d", 0:63)
  centroid <- ee$Array(as.list(model_res$weights))
  scored <- te_fc$map(function(f) {
    pt <- ee$Array(f$toArray(emb_cols))
    f$set("simscore", pt$multiply(centroid)$reduce(ee$Reducer$sum(), list(0L))$get(list(0L)))
  })
  info <- retry_curl_download(scored$select(list(id, "simscore"))$getInfo())
  setNames(vapply(info$features, function(f) as.numeric(f$properties$simscore), numeric(1)),
           vapply(info$features, function(f) as.character(as.integer(f$properties[[id]])), character(1)))
}

#' Sweep methods and hyperparameters and select the best by spatial cross-validation
#'
#' Trains every method x hyperparameter configuration under spatial CV, pools the
#' out-of-fold predictions, and ranks configurations by AUC (or TSS). Works on
#' presence/absence data (the `present` column must contain both 0 and 1). Embeddings
#' are extracted once and reused across the whole grid.
#'
#' @param data Data frame with `longitude`, `latitude`, `present` (0/1), optional `year`.
#' @param methods Methods to tune (subset of svm, rf, gbt, maxent, knn, cart, mindist).
#'   smileNaiveBayes and the presence-only `similarity` reducer are excluded.
#' @param grids Optional named list of per-method configuration lists; default grids
#'   are used for any method not supplied (see `default_tuning_grids`).
#' @param cv_folds,cv_method,block_size Spatial CV settings (blockCV by default).
#' @param scale,aoi_year Embedding sampling scale and year.
#' @param metric Ranking metric, "auc" or "tss".
#' @param materialize If TRUE (default), sample the embeddings ONCE to a temporary GEE
#'   asset and reuse them across the whole grid. Otherwise the embeddings are re-sampled
#'   inside every configuration's `getInfo` (much slower for multi-config sweeps). The
#'   asset is removed when the call finishes.
#' @param gee_project,python_path GEE setup.
#' @return A list with `results` (ranked data frame of method/config/auc/tss, including a
#'   performance-weighted `ensemble` of the best config per method), `best_per_method`,
#'   `best` (overall, which may be the ensemble), and `oof` (pooled OOF scores per config).
#' @export
tune_models <- function(data, methods = c("svm", "rf", "gbt", "maxent", "knn", "cart", "mindist", "similarity"), grids = NULL,
                        cv_folds = 5L, cv_method = "block", block_size = NULL,
                        scale = 10, aoi_year = NULL, metric = c("auc", "tss"),
                        materialize = TRUE,
                        gee_project = NULL, python_path = NULL) {
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)
  ee     <- reticulate::import("ee")
  metric <- match.arg(metric)
  cv_folds <- as.integer(cv_folds)
  if (is.null(aoi_year)) aoi_year <- 2023

  if (!all(c("longitude", "latitude", "present") %in% names(data)))
    stop("`data` must have columns: longitude, latitude, present (0/1).")
  if (length(unique(stats::na.omit(data$present))) < 2)
    stop("`data$present` must contain both presences (1) and absences (0).")
  if (!"year" %in% names(data)) data$year <- aoi_year
  if (is.null(grids)) grids <- default_tuning_grids(methods)

  emb_cols <- sprintf("A%02d", 0:63)
  data$fold <- assign_cv_folds(data, cv_folds, method = cv_method, block_size = block_size)
  data$rid  <- seq_len(nrow(data))

  props <- c("year", "present", "fold", "rid")
  years <- as.list(sort(unique(as.integer(data$year))))
  pts   <- data[, c("longitude", "latitude", props)]

  # Materialise embeddings ONCE so each configuration reads pre-computed features
  # instead of re-running the 64-band sampleRegions (the sweep's real bottleneck).
  assets_to_clean <- character(0)
  on.exit(for (a in assets_to_clean) ee_delete_asset_quietly(a), add = TRUE)
  fc <- NULL
  if (materialize) {
    sdm_section("Materialising embeddings to a GEE asset (sampled once, reused across the grid)")
    res <- tryCatch(
      sample_embeddings_to_asset(pts, scale, properties = props, years = years, project = gee_project),
      error = function(e) { sdm_warn(sprintf("Materialisation failed (%s); re-sampling per configuration.",
                                             conditionMessage(e))); NULL })
    if (!is.null(res)) { fc <- res$fc; assets_to_clean <- res$asset_ids }
  }
  if (is.null(fc)) {
    sdm_section("Sampling embeddings (lazy; re-sampled per configuration)")
    fc <- get_embeddings_at_fc(upload_points_to_gee(pts), scale, properties = props, years = years)
  }
  truth <- setNames(data$present, as.character(data$rid))

  n_cfg <- sum(vapply(methods, function(m) length(grids[[m]]), integer(1)))
  sdm_section(sprintf("Tuning %d configurations across %d method(s), %d-fold %s CV",
                      n_cfg, length(methods), cv_folds,
                      if (cv_method == "block") "blockCV" else "k-means"))
  pb   <- sdm_progress_start("Tuning", total = n_cfg)
  rows <- list(); oof_store <- list()

  for (m in methods) for (cfg in grids[[m]]) {
    oof <- setNames(rep(NA_real_, nrow(data)), as.character(data$rid))
    ok  <- TRUE
    for (k in seq_len(cv_folds)) {
      tr <- fc$filter(ee$Filter$neq("fold", k))$sort("present", FALSE)  # presence-first (C_SVC invert)
      te <- fc$filter(ee$Filter$eq("fold", k))
      sc <- tryCatch({
        mr <- train_gee_model(tr, m, params = cfg, class_property = "present")
        if (isTRUE(mr$is_classifier)) score_test_fc(mr, te) else score_similarity_fc(mr, te)
      }, error = function(e) NULL)
      if (is.null(sc)) { ok <- FALSE; break }
      oof[names(sc)] <- sc
    }
    lbl <- .cfg_label(cfg)
    if (ok && sum(!is.na(oof)) > 0) {
      mm <- calculate_classifier_metrics(oof[truth == 1 & !is.na(oof)], oof[truth == 0 & !is.na(oof)])
      rows[[length(rows) + 1]] <- data.frame(method = m, config = lbl,
        auc = mm$auc_roc, tss = mm$tss, n = sum(!is.na(oof)), stringsAsFactors = FALSE)
      oof_store[[paste(m, lbl)]] <- oof
    } else {
      rows[[length(rows) + 1]] <- data.frame(method = m, config = lbl,
        auc = NA_real_, tss = NA_real_, n = NA_integer_, stringsAsFactors = FALSE)
    }
    pb <- sdm_progress_update(pb)
  }
  sdm_progress_done(pb)

  results <- do.call(rbind, rows)
  keyv    <- if (metric == "auc") results$auc else results$tss
  results <- results[order(-keyv), ]
  best_per_method <- do.call(rbind, lapply(split(results, results$method), function(df) {
    kk <- if (metric == "auc") df$auc else df$tss
    if (all(is.na(kk))) df[1, ] else df[which.max(kk), ]
  }))

  # Performance-WEIGHTED ensemble of the best configuration per method: rank-normalise
  # each method's OOF and combine with weights proportional to its skill above chance
  # (metric - baseline), so weak/near-random methods contribute little. Added as a
  # candidate so it competes with the single configs for "best".
  rids   <- as.character(data$rid)
  bpm_ok <- best_per_method[!is.na(best_per_method$auc), , drop = FALSE]
  if (nrow(bpm_ok) >= 2) {
    keys <- paste(bpm_ok$method, bpm_ok$config)
    rkn  <- function(x) { r <- rank(x, na.last = "keep"); (r - min(r, na.rm = TRUE)) / (max(r, na.rm = TRUE) - min(r, na.rm = TRUE)) }
    R    <- vapply(keys, function(k) rkn(oof_store[[k]][rids]), numeric(length(rids)))
    base <- if (metric == "auc") 0.5 else 0
    w    <- pmax((if (metric == "auc") bpm_ok$auc else bpm_ok$tss) - base, 0)
    if (sum(w) > 0) {
      ens <- setNames(apply(R, 1, function(row) {
        ok <- !is.na(row); if (!any(ok)) NA_real_ else sum(row[ok] * w[ok]) / sum(w[ok]) }), rids)
      tr_e <- truth[rids]; ok2 <- !is.na(ens) & !is.na(tr_e)
      em   <- calculate_classifier_metrics(ens[ok2 & tr_e == 1], ens[ok2 & tr_e == 0])
      results <- rbind(results, data.frame(method = "ensemble", config = "weighted(best-per-method)",
        auc = em$auc_roc, tss = em$tss, n = sum(ok2), stringsAsFactors = FALSE))
      oof_store[["ensemble weighted(best-per-method)"]] <- ens
      keyv    <- if (metric == "auc") results$auc else results$tss
      results <- results[order(-keyv), ]
    }
  }
  best <- results[1, ]

  sdm_section("Tuning complete")
  sdm_done(sprintf("Best: %s [%s]  AUC=%.3f  TSS=%.3f", best$method, best$config, best$auc, best$tss))
  list(results = results, best_per_method = best_per_method, best = best, oof = oof_store)
}
