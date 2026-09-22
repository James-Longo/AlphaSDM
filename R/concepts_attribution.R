# Concept attribution: rank labeled SAE concepts by how much a fitted SDM
# relies on them. Three methods: selection (activation at presences vs
# background), projection (cosine against the linear direction), ablation
# (suitability drop when a concept is removed from the embedding image).

#' Sample concept activations at a set of coordinates
#'
#' Groups the points by year, samples the per-year activation image, and reads
#' the table back. A synchronous read that Earth Engine refuses is re-routed
#' through a batch export, never shrunk.
#'
#' @param df Data frame with longitude, latitude, year.
#' @param sae SAE weights object.
#' @param scale Sampling scale in metres.
#' @param project Earth Engine project id or NULL.
#' @return n x m activation matrix; rows follow `df` (rows whose pixel is
#'   masked, e.g. open water, come back NA and are reported).
#' @noRd
sample_concept_activations <- function(df, sae, scale, project = NULL) {
  ee <- reticulate::import("ee")
  bands <- concept_band_names(sae$m)
  df$row_idx <- seq_len(nrow(df)) - 1L

  fc <- upload_points_to_gee(df[, c("longitude", "latitude", "year", "row_idx")])
  sampled_fcs <- list()
  for (yr in sort(unique(as.integer(df$year)))) {
    act_img <- ee_concept_activations(get_embedding_image(yr), sae)
    yr_fc <- fc$filter(ee$Filter$eq("year", as.integer(yr)))
    sampled_fcs <- c(sampled_fcs, list(
      act_img$sampleRegions(collection = yr_fc, properties = list("row_idx"),
                            scale = scale, geometries = TRUE, tileScale = 16L)))
  }
  sampled <- ee$FeatureCollection(sampled_fcs)$flatten()

  feats <- tryCatch(read_fc_paged(sampled)$features, error = function(e) e)
  if (inherits(feats, "error")) {
    if (!is_gee_timeout(feats)) stop(feats)
    sdm_warn("Activation sampling hit a GEE compute limit; re-routing through a batch export.",
             indent = 1L)
    feats <- ee_table_to_info_async(sampled, project)$features
  }

  out <- matrix(NA_real_, nrow = nrow(df), ncol = sae$m, dimnames = list(NULL, bands))
  for (f in feats) {
    ridx <- suppressWarnings(as.integer(f$properties[["row_idx"]]))
    if (length(ridx) != 1L || is.na(ridx) || ridx < 0L || ridx >= nrow(df)) next
    out[ridx + 1L, ] <- as.numeric(unlist(f$properties[bands]))
  }
  out
}

#' Selection scores: concept activation at presences vs background
#'
#' The score is the ratio of mean activation at presence pixels to mean
#' activation at background pixels, a selection ratio in concept space. The
#' bootstrap resamples rows of both matrices independently. A concept silent in
#' the background but active at presences gives Inf, and one silent in both
#' gives NaN; both are reported as they are rather than smoothed away.
#'
#' @param pres_mat,bg_mat Activation matrices from [sample_concept_activations()].
#' @param boot_reps Bootstrap replicates for the CI.
#' @param conf Confidence level for the bootstrap interval.
#' @return Data frame with per-concept means, ratio and CI.
#' @noRd
selection_scores <- function(pres_mat, bg_mat, boot_reps = 1000L, conf = 0.95) {
  pres_mat <- pres_mat[stats::complete.cases(pres_mat), , drop = FALSE]
  bg_mat   <- bg_mat[stats::complete.cases(bg_mat), , drop = FALSE]
  if (nrow(pres_mat) == 0L || nrow(bg_mat) == 0L)
    stop("No unmasked presence or background activations to score.", call. = FALSE)

  pres_mean <- colMeans(pres_mat)
  bg_mean   <- colMeans(bg_mat)
  ratio     <- pres_mean / bg_mean

  # Bootstrap column means as crossprod with a resample-count vector; this
  # keeps each replicate at one matrix product instead of a full resample.
  boot_ratio <- matrix(NA_real_, nrow = boot_reps, ncol = ncol(pres_mat))
  np <- nrow(pres_mat); nb <- nrow(bg_mat)
  for (b in seq_len(boot_reps)) {
    wp <- tabulate(sample.int(np, np, replace = TRUE), nbins = np)
    wb <- tabulate(sample.int(nb, nb, replace = TRUE), nbins = nb)
    boot_ratio[b, ] <- (crossprod(pres_mat, wp) / np) / (crossprod(bg_mat, wb) / nb)
  }
  # A 0/0 replicate means the concept fired nowhere in that resample: zero
  # evidence of enrichment, counted as 0. Dropping it instead would let a
  # concept supported by a single pixel keep an infinite lower bound.
  boot_ratio[is.nan(boot_ratio)] <- 0
  alpha <- (1 - conf) / 2
  ci <- apply(boot_ratio, 2L, stats::quantile,
              probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)

  data.frame(concept_id = colnames(pres_mat),
             pres_mean = pres_mean, bg_mean = bg_mean,
             pres_active = colMeans(pres_mat > 0), bg_active = colMeans(bg_mat > 0),
             selection_ratio = ratio,
             selection_lo = ci[1L, ], selection_hi = ci[2L, ],
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Projection scores: cosine between a linear model direction and each decoder column
#' @param beta 64-d coefficient vector in the shared embedding value convention.
#' @noRd
projection_scores <- function(beta, sae) {
  beta <- as.numeric(beta)
  if (length(beta) != 64L) stop("beta must be a 64-d vector.", call. = FALSE)
  num <- as.numeric(sae$W_dec %*% beta)
  den <- sqrt(rowSums(sae$W_dec^2)) * sqrt(sum(beta^2))
  data.frame(concept_id = concept_band_names(sae$m),
             projection_cos = num / den,
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Reduce an image's per-band mean over a region, escalating on refusal
#' @noRd
ee_region_means <- function(img, region, scale, project = NULL) {
  ee <- reticulate::import("ee")
  # 1e13 is the largest maxPixels Earth Engine accepts; it exists so the cap
  # never binds and the server's own limits decide.
  dict <- img$reduceRegion(reducer = ee$Reducer$mean(), geometry = region,
                           scale = scale, maxPixels = 1e13, tileScale = 16L)
  res <- tryCatch(retry_curl_download(dict$getInfo()), error = function(e) e)
  if (inherits(res, "error")) {
    if (!is_gee_timeout(res)) stop(res)
    sdm_warn("Region reduction hit a GEE compute limit; re-routing through a batch export.",
             indent = 1L)
    feats <- ee_table_to_info_async(
      ee$FeatureCollection(list(ee$Feature(NULL, dict))), project)$features
    res <- feats[[1L]]$properties
  }
  vapply(res, function(v) if (is.null(v)) NA_real_ else as.numeric(v), numeric(1))
}

#' Ablation scores: mean suitability drop when a concept is removed
#'
#' For concept j the embedding image is perturbed to x - a_j d_j (the
#' subtraction form, so reconstruction error never enters the delta), the SDM
#' prediction is re-run on the perturbed image, and the mean drop over the
#' region is reduced in one request across all requested concepts.
#'
#' @param model_res A fitted model entry from the sdm object's `models`.
#' @param concept_ids Concept band names to ablate.
#' @return list(scores = data.frame, delta_maps = named list of ee.Images).
#' @noRd
ablation_scores <- function(model_res, sae, concept_ids, year, region, scale,
                            project = NULL, maps = FALSE) {
  ee <- reticulate::import("ee")
  emb_cols <- EMB_BANDS
  bands <- concept_band_names(sae$m)

  emb <- get_embedding_image(year)
  act <- ee_concept_activations(emb, sae, select_concepts = concept_ids)
  base_pred <- predict_gee_map(model_res, emb)

  delta_maps <- list()
  delta_img <- NULL
  for (cid in concept_ids) {
    j <- match(cid, bands)
    dk <- ee$Image$constant(as.list(unname(sae$W_dec[j, ])))$rename(emb_cols)
    pert <- emb$subtract(dk$multiply(act$select(cid)))
    delta <- base_pred$subtract(predict_gee_map(model_res, pert))$rename(cid)
    if (maps) delta_maps[[cid]] <- delta
    delta_img <- if (is.null(delta_img)) delta else delta_img$addBands(delta)
  }

  means <- ee_region_means(delta_img, region, scale, project)
  list(scores = data.frame(concept_id = concept_ids,
                           ablation_drop = as.numeric(means[concept_ids]),
                           row.names = NULL, stringsAsFactors = FALSE),
       delta_maps = delta_maps)
}

#' Attribute a fitted SDM to labeled SAE concepts
#'
#' Ranks the concept dictionary (trained once and shipped with the package) by
#' how much the fitted model relies on each concept. Requires the result of
#' [evaluate_models()] or [generate_map()] from the current session: the model
#' handles it carries are live Earth Engine objects and do not survive a
#' restart.
#'
#' Three methods:
#' \describe{
#'   \item{selection}{Mean concept activation at presence pixels over mean
#'     activation at background pixels drawn the same way the SDM drew its
#'     background, with a bootstrap CI. Works for every model.}
#'   \item{projection}{Cosine similarity between the model's 64-d linear
#'     direction and each concept's decoder direction. Only models with a
#'     linear direction qualify (currently `similarity`); runs locally with no
#'     Earth Engine calls.}
#'   \item{ablation}{Re-runs the model's prediction with one concept removed
#'     from the embedding image and reduces the mean suitability drop over the
#'     region. Run for the `top_n` concepts ranked by selection, since each
#'     concept costs a full prediction pass over the region.}
#' }
#'
#' @param sdm Result of [evaluate_models()] or [generate_map()] from this session.
#' @param method Attribution methods to compute: any of `"selection"`,
#'   `"projection"`, `"ablation"`. Defaults to selection.
#' @param top_n Concepts the print method shows, and the number ablation runs on.
#' @param model Which fitted model to attribute, by method name (default: the
#'   first one fitted). Projection additionally requires a linear model.
#' @param year Embedding year for activation sampling and ablation (default:
#'   the year the SDM was fitted on).
#' @param region Region for ablation and background draws (default: the SDM's AOI).
#' @param maps If `TRUE`, attach lazy `ee.Image`s (the concept activation image
#'   and per-concept ablation delta maps) as the `"maps"` attribute.
#' @param boot_reps Bootstrap replicates for the selection CI.
#' @param conf Confidence level for the selection CI.
#' @param sae Optional SAE weights object; default the weights shipped with
#'   the package.
#' @param labels Optional concept label table; default the table shipped with
#'   the package, when present.
#' @return A data frame of class `alphasdm_concepts`: `concept_id`, `label`,
#'   `coherence`, plus score columns for each requested method, ranked by the
#'   first requested method (selection ranks by the bootstrap lower bound of
#'   the ratio, the conservative choice, with mean presence activation as the
#'   tiebreak). Concepts with no label entry are flagged by the print method
#'   so an unlabeled latent is not over-interpreted.
#' @export
derive_concepts <- function(sdm, method = "selection",
                            top_n = 15L, model = NULL, year = NULL, region = NULL,
                            maps = FALSE, boot_reps = 1000L, conf = 0.95,
                            sae = NULL, labels = NULL) {
  methods <- match.arg(method, c("selection", "projection", "ablation"),
                       several.ok = TRUE)
  if (is.null(sdm$models) || is.null(sdm$context)) {
    stop("`sdm` must be the result of evaluate_models() or generate_map() from ",
         "this R session (older results carry no live model handles).", call. = FALSE)
  }
  ctx <- sdm$context
  # projection is pure local matrix math; only the sampled methods touch GEE.
  if (any(c("selection", "ablation") %in% methods))
    ensure_gee_authenticated(project = ctx$gee_project)

  if (is.null(sae)) sae <- load_sae_weights()
  if (is.null(labels)) {
    labels <- tryCatch(concept_labels(), error = function(e) {
      sdm_warn("No concept label table installed; concepts will be reported unlabeled.")
      NULL
    })
  }
  if (is.null(model)) model <- names(sdm$models)[1L]
  if (!model %in% names(sdm$models))
    stop("Model '", model, "' is not in this sdm. Fitted: ",
         paste(names(sdm$models), collapse = ", "), call. = FALSE)
  model_res <- sdm$models[[model]]
  yr  <- if (is.null(year))   ctx$aoi_year else year
  reg <- if (is.null(region)) ctx$aoi_geom else region

  out <- data.frame(concept_id = concept_band_names(sae$m),
                    row.names = NULL, stringsAsFactors = FALSE)
  maps_attr <- list()

  need_selection <- "selection" %in% methods || "ablation" %in% methods
  if (need_selection) {
    sdm_section("Concept attribution: selection")
    pres_df <- ctx$data[ctx$data$present == 1,
                        c("longitude", "latitude", "year"), drop = FALSE]
    n_bg <- nrow(pres_df)
    sdm_info(sprintf("Sampling concept activations at %d presences ...", nrow(pres_df)),
             indent = 1L)
    pres_mat <- sample_concept_activations(pres_df, sae, ctx$scale, ctx$gee_project)

    sdm_info(sprintf("Drawing and sampling %d background points ...", n_bg), indent = 1L)
    # Selection scores compare use against AVAILABILITY, so the background
    # stays a plain random draw; exclusion would change the estimand.
    bg_df <- generate_background_fc_gee(yr, n_bg, reg)$df
    bg_mat <- sample_concept_activations(bg_df, sae, ctx$scale, ctx$gee_project)

    n_na <- sum(!stats::complete.cases(pres_mat)) + sum(!stats::complete.cases(bg_mat))
    if (n_na > 0)
      sdm_warn(sprintf("%d point%s dropped: no satellite coverage.", n_na,
                       if (n_na == 1) "" else "s"), indent = 1L)
    sel <- selection_scores(pres_mat, bg_mat, boot_reps = boot_reps, conf = conf)
    out <- merge(out, sel, by = "concept_id", sort = FALSE)
    if (maps) maps_attr$activations <- ee_concept_activations(get_embedding_image(yr), sae)
  }

  if ("projection" %in% methods) {
    if (isTRUE(model_res$is_classifier)) {
      stop("method = 'projection' needs a model with a linear 64-d direction; '",
           model, "' has none. Fit with methods = 'similarity', or use ",
           "'selection'/'ablation', which work for every model.", call. = FALSE)
    }
    out <- merge(out, projection_scores(model_res$weights, sae),
                 by = "concept_id", sort = FALSE)
  }

  if ("ablation" %in% methods) {
    sdm_section("Concept attribution: ablation")
    ranked <- out$concept_id[order(-out$selection_lo, -out$pres_mean)]
    targets <- ranked[seq_len(min(as.integer(top_n), length(ranked)))]
    sdm_info(sprintf("Ablating the top %d concepts by selection (each is a full %s prediction pass) ...",
                     length(targets), toupper(model)), indent = 1L)
    abl <- ablation_scores(model_res, sae, targets, yr, reg, ctx$scale,
                           ctx$gee_project, maps = maps)
    out <- merge(out, abl$scores, by = "concept_id", all.x = TRUE, sort = FALSE)
    if (maps) maps_attr$ablation_delta <- abl$delta_maps
  }

  if (!is.null(labels)) {
    keep <- intersect(c("concept_id", "label", "coherence"), names(labels))
    out <- merge(out, labels[, keep, drop = FALSE], by = "concept_id",
                 all.x = TRUE, sort = FALSE)
  } else {
    out$label <- NA_character_
    out$coherence <- NA_real_
  }

  # Selection ranks by the bootstrap lower bound, so a ratio inflated by a
  # single supporting pixel cannot outrank a well-supported concept; ties
  # break on mean presence activation.
  ord <- switch(methods[1L],
    selection  = order(-out$selection_lo, -out$pres_mean),
    projection = order(-out$projection_cos),
    ablation   = order(-out$ablation_drop))
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  front <- intersect(c("concept_id", "label", "coherence"), names(out))
  out <- out[, c(front, setdiff(names(out), front)), drop = FALSE]

  class(out) <- c("alphasdm_concepts", "data.frame")
  attr(out, "methods") <- methods
  attr(out, "model") <- model
  attr(out, "top_n") <- as.integer(top_n)
  if (maps) attr(out, "maps") <- maps_attr
  out
}

#' @export
print.alphasdm_concepts <- function(x, n = NULL, ...) {
  if (is.null(n)) n <- attr(x, "top_n")
  if (is.null(n)) n <- nrow(x)
  n <- min(n, nrow(x))
  cat(sprintf("Concept attribution (%s) for model '%s': top %d of %d concepts\n",
              paste(attr(x, "methods"), collapse = " + "), attr(x, "model"),
              n, nrow(x)))
  df <- as.data.frame(x)[seq_len(n), , drop = FALSE]
  unlabeled <- is.na(df$label)
  df$label[unlabeled] <- "(unlabeled latent)"
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(v) signif(v, 3))
  print.data.frame(df, row.names = FALSE)
  if (any(unlabeled)) {
    cat("! Unlabeled concepts have no catalog correlate on record;",
        "do not over-interpret them.\n")
  }
  if (!is.null(attr(x, "maps")))
    cat("Lazy ee.Image maps attached; retrieve with attr(x, \"maps\").\n")
  invisible(x)
}
