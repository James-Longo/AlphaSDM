# Pseudo-absence generation following Barbet-Massin et al. (2012 MEE).
# Two pools serve the two method groups the paper distinguishes:
#   - the FULL pool (maxent, svm, cart, ...): a large plain-random draw,
#     10,000 points by default (their recommendation, endorsing Phillips &
#     Dudik 2008), never subjected to exclusion — random placement is what
#     the paper recommends for these methods;
#   - the BALANCED pool (rf, gbt, knn): as many points as presences, and
#     optionally drawn OUTSIDE both the geographic neighbourhood and the
#     environmental envelope of the presences (their "2-degree-far" + SRE
#     recipes for classification methods), with both thresholds estimated
#     from the presence data rather than fixed:
#       geographic radius   = embedding-autocorrelation range
#       environmental edge  = bias-corrected max Mahalanobis distance of
#                             the presences to their own embedding cloud
#     (a min/max box, the literal SRE, excludes nothing in 64 dims).
# Every pool tops itself up until it reaches its target: coverage losses
# (water, masked pixels) and exclusion losses trigger further draws, and
# the geographic radius relaxes stepwise before the target is ever missed,
# because an under-sized background changes prevalence, which the paper
# shows drives model accuracy directly.

# Computational constants (not ecological choices): the similarity-decay
# distance ladder, pair/bootstrap counts, oversampling factor, top-up
# rounds.
BG_LADDER_M <- 250 * 2^(0:7)   # 250 m .. 32 km
BG_MAX_PAIRS <- 200L
BG_BOOT <- 200L
BG_OVERSAMPLE <- 3L
BG_TOPUP_ROUNDS <- 5L
BG_ENV_CAP <- 2000L            # presence rows used for the envelope
BG_DEFAULT_POOL <- 10000L      # full-pool default (Barbet-Massin 2012)
BG_DRAW_CHUNK <- 5000L         # points per probe request (adaptive)

#' Embedding-autocorrelation exclusion radius.
#'
#' Samples embeddings at random-bearing offsets on a doubling distance
#' ladder around up to BG_MAX_PAIRS presences. The baseline is the
#' presence-to-random-REGIONAL similarity (presence-to-presence pairs are
#' self-similar by habitat selection and would collapse the radius). The
#' radius is the smallest ladder distance whose mean presence-offset
#' cosine similarity is statistically indistinguishable (bootstrap 95% CI)
#' from that baseline.
#' @return list(radius_m, decay, baseline)
#' @noRd
estimate_similarity_range <- function(pres_df, region, aoi_year, scale,
                                      seed = 0L) {
  ee <- reticulate::import("ee")
  set.seed(seed + 1L)
  np <- min(BG_MAX_PAIRS, nrow(pres_df))
  pd <- pres_df[sample.int(nrow(pres_df), np), c("longitude", "latitude")]

  grid <- expand.grid(i = seq_len(np), d = BG_LADDER_M)
  grid$theta <- runif(nrow(grid), 0, 2 * pi)
  grid$latitude <- pd$latitude[grid$i] + grid$d * cos(grid$theta) / 111320
  grid$longitude <- pd$longitude[grid$i] +
    grid$d * sin(grid$theta) / (111320 * cos(pd$latitude[grid$i] * pi / 180))

  up <- rbind(
    data.frame(longitude = pd$longitude, latitude = pd$latitude,
               pair = seq_len(np), dist = 0),
    data.frame(longitude = grid$longitude, latitude = grid$latitude,
               pair = grid$i, dist = grid$d))
  up$year <- as.integer(aoi_year)
  samp <- get_embeddings_at_fc(upload_points_to_gee(up), scale,
                               properties = c("pair", "dist"),
                               years = list(as.integer(aoi_year)))
  rows <- read_fc_paged(samp)$features
  emb_cols <- sprintf("A%02d", 0:63)
  val <- as.data.frame(do.call(rbind, lapply(rows, function(f) {
    p <- f$properties
    c(pair = as.numeric(p$pair), dist = as.numeric(p$dist),
      vapply(emb_cols, function(cn)
        if (is.null(p[[cn]])) NA_real_ else as.numeric(p[[cn]]), numeric(1)))
  })))
  E <- as.matrix(val[, emb_cols])
  E <- E / pmax(sqrt(rowSums(E^2)), 1e-12)

  base_idx <- which(val$dist == 0)
  anchors <- E[base_idx, , drop = FALSE]
  names(base_idx) <- val$pair[base_idx]

  rnd_fc <- ee$FeatureCollection$randomPoints(
    region, as.integer(2L * BG_MAX_PAIRS), as.integer(seed))$
    map(function(f) f$set("year", as.integer(aoi_year)))
  rnd_samp <- get_embeddings_at_fc(rnd_fc, scale, properties = list("year"),
                                   years = list(as.integer(aoi_year)))
  R <- as.matrix(as.data.frame(do.call(rbind,
    lapply(read_fc_paged(rnd_samp)$features, function(f)
      vapply(emb_cols, function(cn) {
        v <- f$properties[[cn]]
        if (is.null(v)) NA_real_ else as.numeric(v)
      }, numeric(1))))))
  R <- R[stats::complete.cases(R), , drop = FALSE]
  R <- R / pmax(sqrt(rowSums(R^2)), 1e-12)
  ri <- sample.int(nrow(R), nrow(anchors), replace = nrow(R) < nrow(anchors))
  baseline <- mean(rowSums(anchors * R[ri, , drop = FALSE]))

  decay <- do.call(rbind, lapply(BG_LADDER_M, function(d) {
    ii <- which(val$dist == d)
    if (!length(ii)) return(NULL)
    a <- anchors[match(as.character(val$pair[ii]), names(base_idx)), ,
                 drop = FALSE]
    ok <- stats::complete.cases(a) & stats::complete.cases(E[ii, , drop = FALSE])
    sims <- rowSums(a[ok, , drop = FALSE] * E[ii[ok], , drop = FALSE])
    if (length(sims) < 10) return(NULL)
    bs <- replicate(BG_BOOT, mean(sample(sims, replace = TRUE)))
    data.frame(dist = d, sim = mean(sims),
               lo = stats::quantile(bs, 0.025), hi = stats::quantile(bs, 0.975))
  }))
  hit <- decay$dist[decay$lo <= baseline]
  radius <- if (length(hit)) min(hit) else max(BG_LADDER_M)
  list(radius_m = radius, decay = decay, baseline = baseline)
}

#' Bias-corrected Mahalanobis envelope of the presence embeddings.
#'
#' Diagonal-covariance Mahalanobis distance to the presence cloud; the
#' exclusion threshold is the bootstrap bias-corrected maximum
#' (2*obs_max - mean(bootstrap maxima)). Few presences produce a larger
#' correction and a wider exclusion zone, replacing the paper's
#' few-vs-many presence-count switch with a continuous rule.
#' @noRd
estimate_embedding_envelope <- function(P, seed = 0L) {
  set.seed(seed + 2L)
  if (nrow(P) > BG_ENV_CAP) P <- P[sample.int(nrow(P), BG_ENV_CAP), ]
  mu <- colMeans(P); sdv <- pmax(apply(P, 2, stats::sd), 1e-9)
  d <- sqrt(rowSums(sweep(sweep(P, 2, mu), 2, sdv, "/")^2))
  obs_max <- max(d)
  boot_max <- replicate(BG_BOOT, max(sample(d, replace = TRUE)))
  list(mu = mu, sd = sdv, threshold = 2 * obs_max - mean(boot_max))
}


#' Generates a background pool on Earth Engine, guaranteed to reach its
#' target size.
#'
#' All strategies share one structure: gather valid-coverage candidate
#' coordinates in bounded chunks (one interactive request over a large
#' pool times out), filter them by the active exclusions, upload the kept
#' coordinates, and sample the 64 bands once. Geographic exclusion halves
#' its radius stepwise before a shortfall is ever accepted; the
#' environmental envelope (the false-absence guard) never relaxes.
#' @param use_geo Draw outside a buffer around the presences; radius from
#'   `radius_m` or, when NULL, the embedding-autocorrelation range.
#' @param use_env Drop candidates inside the presence Mahalanobis
#'   envelope; threshold from `env_threshold` or, when NULL, the
#'   bias-corrected presence maximum.
#' @param seed Integer; varies the draw for replicate background sets.
#' @param est Precomputed list(rng, env) to reuse across replicate draws.
#' @return list(df, fc, n_returned, radius_m, env_threshold,
#'   n_candidates, n_env_excluded, est)
#' @noRd
generate_background_fc_gee <- function(aoi_year, count, region, scale = 10,
                                       presence_df = NULL,
                                       presence_emb = NULL,
                                       use_geo = FALSE, use_env = FALSE,
                                       radius_m = NULL, env_threshold = NULL,
                                       seed = 0L, est = NULL) {
  ee <- reticulate::import("ee")
  count <- as.integer(count)
  emb_cols <- sprintf("A%02d", 0:63)

  rng <- NULL; env <- NULL
  if (use_geo && is.null(radius_m)) {
    if (is.null(est)) {
      if (is.null(presence_df))
        stop("Geographic exclusion needs presence coordinates to estimate ",
             "the radius; pass `radius_m` to set it directly.", call. = FALSE)
      rng <- estimate_similarity_range(presence_df, region, aoi_year, scale,
                                       seed = seed)
    } else rng <- est$rng
    radius_m <- rng$radius_m
    sdm_info(sprintf(
      "Exclusion radius %.0f m (similarity decay; baseline cos %.3f)",
      radius_m, rng$baseline), indent = 1L)
  }
  if (use_env && is.null(env_threshold)) {
    if (is.null(est) || is.null(est$env)) {
      if (is.null(presence_emb))
        stop("Envelope exclusion needs presence embeddings to estimate the ",
             "threshold; pass `env_threshold` to set it directly.",
             call. = FALSE)
      env <- estimate_embedding_envelope(as.matrix(presence_emb), seed = seed)
    } else env <- est$env
    env_threshold <- env$threshold
    sdm_info(sprintf(
      "Envelope: Mahalanobis threshold %.1f (bias-corrected presence max)",
      env_threshold), indent = 1L)
  }
  if (use_env && is.null(env)) {
    # A user-supplied threshold still needs the presence mu/sd frame.
    if (is.null(presence_emb))
      stop("Envelope exclusion needs presence embeddings.", call. = FALSE)
    P <- as.matrix(presence_emb)
    env <- list(mu = colMeans(P), sd = pmax(apply(P, 2, stats::sd), 1e-9),
                threshold = env_threshold)
  }

  bg_region <- region
  if (use_geo) {
    pres_fc <- upload_points_to_gee(
      presence_df[, c("longitude", "latitude")])
    bg_region <- region$difference(pres_fc$geometry()$buffer(radius_m), 100)
  }

  probe_img <- alphaearth_rescale(get_embedding_image(aoi_year))
  probe_1band <- probe_img$select("A00")

  kept <- NULL; n_cand <- 0L; n_excl <- 0L
  per_draw <- BG_DRAW_CHUNK
  max_rounds <- as.integer(ceiling(count / min(per_draw, count)) *
                           BG_TOPUP_ROUNDS * 2L)
  for (r in seq_len(max_rounds)) {
    got <- if (is.null(kept)) 0L else nrow(kept)
    if (got >= count) break
    feats <- tryCatch({
      raw <- ee$FeatureCollection$randomPoints(
        bg_region, as.integer(min(per_draw, count)),
        as.integer(seed + 1000L * (r - 1L)))
      # Candidates carry only the bands the filter needs: one band when no
      # envelope is active, all 64 when it is.
      img_r <- if (use_env) probe_img else probe_1band
      read_fc_paged(img_r$sampleRegions(collection = raw, scale = scale,
                                        geometries = TRUE,
                                        tileScale = 16L))$features
    }, error = function(e) {
      sdm_info(sprintf("Draw failed (%s); adapting",
                       substr(conditionMessage(e), 1, 60)), indent = 1L)
      NULL
    })
    if (is.null(feats)) {
      per_draw <- max(1000L, per_draw %/% 2L)
      next
    }
    need_cols <- if (use_env) emb_cols else "A00"
    cand <- do.call(rbind, lapply(feats, function(f) {
      cc <- f$geometry$coordinates
      c(as.numeric(cc[[1]]), as.numeric(cc[[2]]),
        vapply(need_cols, function(cn) {
          v <- f$properties[[cn]]
          if (is.null(v)) NA_real_ else as.numeric(v)
        }, numeric(1)))
    }))
    if (!is.null(cand)) {
      cand <- as.data.frame(cand)
      names(cand) <- c("longitude", "latitude", need_cols)
      cand <- cand[stats::complete.cases(cand), ]
      n_cand <- n_cand + nrow(cand)
      if (use_env && nrow(cand)) {
        dc <- sqrt(rowSums(sweep(sweep(as.matrix(cand[, emb_cols]), 2,
                                       env$mu), 2, env$sd, "/")^2))
        n_excl <- n_excl + sum(dc <= env$threshold)
        cand <- cand[dc > env$threshold, , drop = FALSE]
      }
      kept <- unique(rbind(kept, cand[, c("longitude", "latitude")]))
    }
    got <- if (is.null(kept)) 0L else nrow(kept)
    if (got >= count) break
    # Geographic exclusion relaxes before a shortfall is accepted; the
    # envelope stays.
    if (use_geo && (r >= 2L || got == 0L) && got < count / 2 &&
        radius_m > min(BG_LADDER_M)) {
      radius_m <- radius_m / 2
      sdm_info(sprintf(
        "Supply short (%d of %d); halving the geographic radius to %.0f m",
        got, count, radius_m), indent = 1L)
      bg_region <- region$difference(
        pres_fc$geometry()$buffer(radius_m), 100)
    } else if (r %% 4L == 0L) {
      sdm_info(sprintf("Background supply %d of %d (round %d)",
                       got, count, r), indent = 1L)
    }
  }
  if (is.null(kept) || nrow(kept) == 0L)
    stop("No background point survived; the AOI appears to be entirely ",
         "presence-like habitat or without coverage. Try strategy = ",
         "\"random\" or a different AOI.", call. = FALSE)
  if (nrow(kept) < count)
    sdm_warn(sprintf(
      "Region supplied only %d of %d background points",
      nrow(kept), count), indent = 1L)
  set.seed(seed + 3L)
  if (nrow(kept) > count) kept <- kept[sample.int(nrow(kept), count), ]
  kept$year <- as.integer(aoi_year); kept$present <- 0L
  out_fc <- get_embeddings_at_fc(upload_points_to_gee(kept), scale,
                                 properties = c("year", "present"),
                                 geometries = TRUE,
                                 years = list(as.integer(aoi_year)))
  list(df = kept, fc = out_fc, n_returned = nrow(kept),
       radius_m = if (use_geo) radius_m else NA_real_,
       env_threshold = if (use_env) env$threshold else NA_real_,
       n_candidates = n_cand, n_env_excluded = n_excl,
       est = list(rng = rng, env = env))
}
