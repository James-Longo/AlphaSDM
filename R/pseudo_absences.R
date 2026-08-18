#' Generate pseudo-absences for presence-only data
#'
#' Presence-only records cannot be modelled directly: every method in
#' AlphaSDM needs absence or background data, and how that background is
#' placed is a modelling decision with real consequences (Barbet-Massin
#' et al. 2012, *Methods in Ecology and Evolution* 3:327-338). This
#' function makes that decision explicit. It draws pseudo-absences inside
#' `aoi` under the strategy you choose, reports every threshold it used,
#' and returns your presences and the new absences as one data frame
#' ready for [evaluate_models()] or [generate_map()].
#'
#' Strategies, following Barbet-Massin et al. (2012):
#' \describe{
#'   \item{`"random"`}{Uniform over the AOI. Their recommendation for
#'     regression-style methods and MaxEnt (with `n = 10000`). MaxEnt is
#'     not in the default ensemble for exactly this reason: give it its
#'     own random set and run `methods = "maxent"` separately.}
#'   \item{`"disk"`}{Their "2-degree-far": only beyond a distance from
#'     every presence. `radius_m = NULL` estimates the distance at which
#'     embedding similarity to the presences decays to the regional
#'     baseline, and reports it; pass a number to choose it yourself.}
#'   \item{`"envelope"`}{Their SRE, in embedding space: only outside the
#'     presence environmental envelope, measured as Mahalanobis distance
#'     to the presence cloud. `env_threshold = NULL` uses the
#'     bias-corrected maximum distance among the presences themselves.}
#'   \item{`"combined"`}{Both exclusions at once. Their recommendation
#'     for classification and machine-learning methods (rf, gbt) with
#'     `n` near the number of presences; validated here on Bicknell's
#'     Thrush (real-absence AUC) and *Prunus africana* (Boyce index).}
#' }
#'
#' Supply is guaranteed: every strategy redraws until `n` points with
#' satellite coverage survive the active exclusions; `"disk"` and
#' `"combined"` halve the radius stepwise rather than come up short (the
#' envelope never relaxes, since points inside it are the likely false
#' absences the strategy exists to avoid).
#'
#' @param data Formatted presence records from [format_data()] (standard
#'   `longitude`, `latitude`, `year`, `present` columns, all presences).
#'   The workflow is format first, then add absences:
#'   \preformatted{
#'   pres <- format_data(obs, coords = c("lon", "lat"), year = "yr")
#'   data <- generate_pseudo_absences(pres, aoi = ..., strategy = ...)
#'   evaluate_models(data)
#'   }
#' @param aoi Where absences may be placed: an `ee.Geometry`, a
#'   `list(lon, lat, radius)`, a path to a vector file, or the string
#'   `"bbox"` to use the presence bounding box (an explicit choice, not a
#'   silent default — a bounding box is rarely the right availability
#'   frame for clustered records).
#' @param strategy One of `"random"`, `"disk"`, `"envelope"`,
#'   `"combined"`. No default: this is the modelling decision.
#' @param n Number of pseudo-absences (default 10000, Barbet-Massin et
#'   al. 2012; use about the presence count for `"combined"` feeding
#'   tree methods).
#' @param radius_m Disk radius in metres; NULL estimates it from the
#'   embedding-autocorrelation range and reports it.
#' @param env_threshold Mahalanobis envelope threshold; NULL uses the
#'   bias-corrected presence maximum and reports it.
#' @param aoi_year Embedding year for placement checks (default: latest
#'   Alpha Earth year).
#' @param scale Sampling scale in metres (default 10).
#' @param seed Integer seed for the draws.
#' @param gee_project Optional Earth Engine cloud project.
#' @return A data frame with `longitude`, `latitude`, `year`, `present`
#'   (your presences as 1, pseudo-absences as 0), ready for
#'   [evaluate_models()] or [generate_map()] directly, carrying the
#'   settings used in `attr(, "pa_settings")`.
#' @export
generate_pseudo_absences <- function(data, aoi, strategy,
                                     n = 10000L,
                                     radius_m = NULL, env_threshold = NULL,
                                     aoi_year = NULL, scale = 10,
                                     seed = 0L, gee_project = NULL) {
  if (missing(strategy))
    stop("Choose a pseudo-absence `strategy`: \"random\", \"disk\", ",
         "\"envelope\" or \"combined\". This is a modelling decision ",
         "AlphaSDM will not make for you; see ?generate_pseudo_absences ",
         "for the recipe per method (Barbet-Massin et al. 2012).",
         call. = FALSE)
  strategy <- match.arg(strategy,
                        c("random", "disk", "envelope", "combined"))
  if (missing(aoi))
    stop("Choose an `aoi` for the pseudo-absences: an ee.Geometry, ",
         "list(lon, lat, radius), a vector file path, or \"bbox\" to use ",
         "the presence bounding box explicitly.", call. = FALSE)
  need <- c("longitude", "latitude", "year", "present")
  if (!is.data.frame(data) || !all(need %in% names(data)))
    stop("`data` must come from format_data() (columns ",
         paste(need, collapse = ", "), "). Format first:\n  pres <- ",
         "format_data(obs, coords = c(...), year = ...)\n  data <- ",
         "generate_pseudo_absences(pres, aoi = ..., strategy = ...)",
         call. = FALSE)
  if (any(data$present == 0))
    stop("`data` already contains absence rows; pseudo-absences are for ",
         "presence-only data.", call. = FALSE)
  if (!is.null(gee_project)) gee_project <- as.character(gee_project)
  ensure_gee_authenticated(project = gee_project)
  ee <- reticulate::import("ee")

  pres <- data
  yrs <- alphaearth_year_range()
  if (is.null(aoi_year)) aoi_year <- yrs[2]

  aoi_geom <- if (identical(aoi, "bbox")) {
    ee$Geometry$Rectangle(c(min(pres$longitude), min(pres$latitude),
                            max(pres$longitude), max(pres$latitude)))
  } else resolve_aoi(aoi, ee)

  use_geo <- strategy %in% c("disk", "combined")
  use_env <- strategy %in% c("envelope", "combined")

  pres_emb <- NULL
  if (use_env) {
    sdm_info("Sampling presence embeddings for the envelope ...",
             indent = 1L)
    pfc <- upload_points_to_gee(
      transform(pres[, c("longitude", "latitude", "year")],
                present = 1L))
    ps <- get_embeddings_at_fc(pfc, scale,
                               properties = c("year", "present"),
                               years = as.list(unique(as.integer(
                                 c(pres$year, aoi_year)))))
    pe <- read_fc_paged(ps$limit(2000L))$features
    emb_cols <- sprintf("A%02d", 0:63)
    pres_emb <- do.call(rbind, lapply(pe, function(f)
      vapply(emb_cols, function(cn) {
        v <- f$properties[[cn]]
        if (is.null(v)) NA_real_ else as.numeric(v)
      }, numeric(1))))
    pres_emb <- pres_emb[stats::complete.cases(pres_emb), , drop = FALSE]
  }

  sdm_section(sprintf("Generating %d pseudo-absences (%s)", n, strategy))
  bg <- generate_background_fc_gee(
    aoi_year, n, aoi_geom, scale = scale,
    presence_df = pres, presence_emb = pres_emb,
    use_geo = use_geo, use_env = use_env,
    radius_m = radius_m, env_threshold = env_threshold, seed = seed)

  sdm_done(sprintf(
    "%d pseudo-absences (%s%s%s)", bg$n_returned, strategy,
    if (use_geo) sprintf("; radius %.0f m", bg$radius_m) else "",
    if (use_env) sprintf("; Mahalanobis > %.1f", bg$env_threshold) else ""))

  abs_df <- bg$df[, c("longitude", "latitude", "year")]
  # A random draw can land on a presence pixel; drop exact collisions so
  # the same location never appears with both labels.
  key <- function(d) paste(round(d$longitude, 6), round(d$latitude, 6))
  dup <- key(abs_df) %in% key(pres)
  if (any(dup)) {
    sdm_warn(sprintf("%d pseudo-absence%s dropped: same location as a presence",
                     sum(dup), if (sum(dup) == 1) "" else "s"), indent = 1L)
    abs_df <- abs_df[!dup, , drop = FALSE]
  }
  out <- rbind(
    transform(pres[, c("longitude", "latitude", "year")], present = 1L),
    transform(abs_df, present = 0L))
  rownames(out) <- NULL
  attr(out, "pa_settings") <- list(
    strategy = strategy, n = bg$n_returned, radius_m = bg$radius_m,
    env_threshold = bg$env_threshold, aoi_year = aoi_year, seed = seed,
    n_candidates = bg$n_candidates, n_env_excluded = bg$n_env_excluded)
  out
}
