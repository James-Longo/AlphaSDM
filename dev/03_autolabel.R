# Phase 4: auto-label SAE latents against GEE catalog layers. Run once per SAE
# version, by hand. Produces inst/extdata/concept_labels_v1.csv and a
# coherence comparison against the NMF and k-means baselines.
#
# Per concept: Spearman correlation with continuous layers, per-class mean
# activation and top-activation enrichment for categorical layers, a coherence
# score (the strongest normalized association), a rule-based label from the
# top correlates, and the top-activating locations for chip verification.
#
# Verify every catalog asset id below at run time; ids drift between releases.

library(AlphaSDM)

# ---- Config ------------------------------------------------------------------
SAE_PATH     <- "dev/artifacts/sae_candidates/sae_topk_m512_k16.rds"  # candidate under test
SAMPLE_DIR   <- "~/alphasdm_sae_sample_v1"   # Phase 1 CSVs (locations reused)
OUT_CSV      <- "dev/artifacts/concept_labels_candidate.csv"
N_LABEL      <- 200000L                       # labeling subsample of Phase 1 points
N_CHIPS      <- 20L                           # example locations kept per concept
TOP_SHARE    <- 0.01                          # "top-activating" = top 1% of pixels
YEAR         <- 2023
SCALE        <- 10
CHUNK        <- 5000L                         # points per sampleRegions request

AlphaSDM:::ensure_gee_authenticated()
ee <- reticulate::import("ee")
sae <- AlphaSDM:::validate_sae(readRDS(SAE_PATH))
bands <- AlphaSDM:::concept_band_names(sae$m)

# ---- Candidate meaning layers ------------------------------------------------
# name = band name in the sampled table; type governs the statistic.
terrain <- ee$Terrain$products(ee$Image("USGS/SRTMGL1_003"))
tc <- ee$ImageCollection("IDAHO_EPSCOR/TERRACLIMATE")$
  filterDate(sprintf("%d-01-01", YEAR), sprintf("%d-01-01", YEAR + 1))$mean()
dw <- ee$ImageCollection("GOOGLE/DYNAMICWORLD/V1")$
  filterDate(sprintf("%d-01-01", YEAR), sprintf("%d-01-01", YEAR + 1))$mean()
viirs <- ee$ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG")$
  filterDate(sprintf("%d-01-01", YEAR), sprintf("%d-01-01", YEAR + 1))$
  select("avg_rad")$median()
ndvi <- ee$ImageCollection("MODIS/061/MOD13Q1")$
  filterDate(sprintf("%d-01-01", YEAR), sprintf("%d-01-01", YEAR + 1))$select("NDVI")
hansen <- ee$Image("UMD/hansen/global_forest_change_2024_v1_12")   # verify version

layers <- list(
  # Climate
  bio01 = list(img = ee$Image("WORLDCLIM/V1/BIO")$select("bio01"), type = "continuous", desc = "annual mean temperature"),
  bio12 = list(img = ee$Image("WORLDCLIM/V1/BIO")$select("bio12"), type = "continuous", desc = "annual precipitation"),
  aet   = list(img = tc$select("aet"), type = "continuous", desc = "actual evapotranspiration"),
  soil_moisture = list(img = tc$select("soil"), type = "continuous", desc = "soil moisture"),
  # Structure
  canopy_height = list(img = ee$Image("users/nlang/ETH_GlobalCanopyHeight_2020_10m_v1"),  # verify id
                       type = "continuous", desc = "canopy height"),
  # Land cover
  worldcover = list(img = ee$ImageCollection("ESA/WorldCover/v200")$first()$select("Map"),
                    type = "categorical", desc = "WorldCover class"),
  dw_trees = list(img = dw$select("trees"), type = "continuous", desc = "Dynamic World tree probability"),
  dw_water = list(img = dw$select("water"), type = "continuous", desc = "Dynamic World water probability"),
  dw_crops = list(img = dw$select("crops"), type = "continuous", desc = "Dynamic World crop probability"),
  dw_built = list(img = dw$select("built"), type = "continuous", desc = "Dynamic World built probability"),
  # Terrain
  elevation = list(img = terrain$select("elevation"), type = "continuous", desc = "elevation"),
  slope     = list(img = terrain$select("slope"), type = "continuous", desc = "slope"),
  # Disturbance
  forest_loss = list(img = hansen$select("lossyear")$gt(0), type = "continuous", desc = "forest loss 2001+"),
  # Human
  nightlights = list(img = viirs, type = "continuous", desc = "VIIRS nightlights"),
  # Phenology / greenness
  ndvi_peak = list(img = ndvi$max()$multiply(0.0001), type = "continuous", desc = "peak NDVI"),
  ndvi_amp  = list(img = ndvi$max()$subtract(ndvi$min())$multiply(0.0001), type = "continuous", desc = "NDVI amplitude")
)

WC_NAMES <- c(`10` = "tree cover", `20` = "shrubland", `30` = "grassland",
              `40` = "cropland", `50` = "built-up", `60` = "bare/sparse",
              `70` = "snow and ice", `80` = "water", `90` = "wetland",
              `95` = "mangrove", `100` = "moss and lichen")

# ---- Labeling sample: reuse Phase 1 locations --------------------------------
files <- list.files(path.expand(SAMPLE_DIR), pattern = "\\.csv$", full.names = TRUE)
pts <- do.call(rbind, lapply(files, function(f)
  utils::read.csv(f)[, c("longitude", "latitude")]))
set.seed(42)
pts <- pts[sample.int(nrow(pts), min(N_LABEL, nrow(pts))), ]
pts$year <- YEAR

# ---- Sample activations and all layers in one pipeline -----------------------
layer_imgs <- Reduce(function(a, b) a$addBands(b),
                     lapply(names(layers), function(nm) layers[[nm]]$img$rename(nm)))
act_img <- AlphaSDM:::ee_concept_activations(
  AlphaSDM:::get_embedding_image(YEAR), sae)
stack <- act_img$addBands(layer_imgs)$addBands(ee$Image$pixelLonLat())

rows <- list()
for (s in seq(1L, nrow(pts), by = CHUNK)) {
  chunk <- pts[s:min(s + CHUNK - 1L, nrow(pts)), ]
  fc <- AlphaSDM:::upload_points_to_gee(chunk)
  sampled <- stack$sampleRegions(collection = fc, properties = list(),
                                 scale = SCALE, geometries = FALSE, tileScale = 16L)
  feats <- tryCatch(AlphaSDM:::read_fc_paged(sampled)$features, error = function(e) e)
  if (inherits(feats, "error")) {
    if (!AlphaSDM:::is_gee_timeout(feats)) stop(feats)
    feats <- AlphaSDM:::ee_table_to_info_async(sampled)$features
  }
  rows <- c(rows, lapply(feats, function(f) as.data.frame(f$properties)))
  message(sprintf("%d / %d points sampled", min(s + CHUNK - 1L, nrow(pts)), nrow(pts)))
}
tab <- do.call(rbind, rows)
saveRDS(tab, "dev/artifacts/autolabel_sample.rds")

# ---- Score each concept ------------------------------------------------------
cont_layers <- names(layers)[vapply(layers, function(l) l$type == "continuous", logical(1))]
label_rows <- list()
for (cid in bands) {
  a <- tab[[cid]]
  if (all(is.na(a)) || stats::sd(a, na.rm = TRUE) == 0) {
    label_rows[[cid]] <- data.frame(concept_id = cid, label = NA, coherence = 0,
                                    dead = TRUE, correlates = "", example_locations = "")
    next
  }
  # Continuous: Spearman correlation.
  cors <- vapply(cont_layers, function(nm)
    suppressWarnings(stats::cor(a, tab[[nm]], method = "spearman",
                                use = "pairwise.complete.obs")), numeric(1))
  # Categorical (WorldCover): enrichment of the top-activating pixels.
  top_cut <- stats::quantile(a, 1 - TOP_SHARE, na.rm = TRUE)
  top <- !is.na(a) & a >= top_cut & a > 0
  cls <- as.character(tab$worldcover)
  enrich <- vapply(unique(stats::na.omit(cls)), function(cl) {
    (mean(cls[top] == cl, na.rm = TRUE)) / max(mean(cls == cl, na.rm = TRUE), 1e-12)
  }, numeric(1))
  best_cls <- names(enrich)[which.max(enrich)]

  # Coherence: strongest normalized association across all layers. Enrichment
  # is mapped to [0, 1] as 1 - 1/e so a 1x (chance) enrichment scores 0.
  best_cor <- max(abs(cors), na.rm = TRUE)
  best_enr <- max(enrich, na.rm = TRUE)
  coherence <- max(best_cor, 1 - 1 / max(best_enr, 1))

  ord <- order(-abs(cors))
  top_correlates <- paste(sprintf("%s=%.2f", cont_layers[ord][1:3], cors[ord][1:3]),
                          collapse = "; ")
  correlates <- sprintf("%s | top-class %s (%s) %.1fx", top_correlates,
                        best_cls, WC_NAMES[best_cls], best_enr)

  # Rule-based label: dominant class plus the strongest continuous qualifiers.
  quals <- cont_layers[ord][1:2]
  dirs  <- ifelse(cors[ord][1:2] > 0, "high", "low")
  label <- sprintf("%s, %s %s + %s %s", WC_NAMES[best_cls],
                   dirs[1], layers[[quals[1]]]$desc, dirs[2], layers[[quals[2]]]$desc)

  ord_a <- order(-a)
  ex <- tab[ord_a[seq_len(min(N_CHIPS, sum(!is.na(a))))], c("longitude", "latitude")]
  label_rows[[cid]] <- data.frame(
    concept_id = cid, label = label, coherence = coherence, dead = FALSE,
    correlates = correlates,
    example_locations = paste(sprintf("%.5f,%.5f", ex$longitude, ex$latitude),
                              collapse = ";"))
}
labels <- do.call(rbind, label_rows)
utils::write.csv(labels, OUT_CSV, row.names = FALSE)

message(sprintf("Mean coherence: %.3f (dead: %d/%d)",
                mean(labels$coherence[!labels$dead]), sum(labels$dead), sae$m))
message("Run this script over each SAE candidate AND the NMF / k-means baseline")
message("dictionaries (wrap a baseline as an SAE via AlphaSDM:::new_sae with")
message("W_enc = W_dec = the dictionary, b = 0, k = m). Compare mean coherence;")
message("if the SAE does not beat the baselines, ship NMF and say so in the vignette.")
message(sprintf("Freeze the winner: copy its rds to inst/extdata/sae_weights_v1.rds "))
message(sprintf("and this CSV to inst/extdata/concept_labels_v1.csv."))
