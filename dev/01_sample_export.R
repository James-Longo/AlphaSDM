# Phase 1: export the SAE training sample from Earth Engine. Run once, by hand.
#
# Draws a globally stratified pixel sample of Alpha Earth embeddings: roughly
# equal draws per ESA WorldCover class per continent, so rare classes such as
# mangrove and wetland are represented rather than swamped by cropland. One
# batch export task per continent, to Google Drive as CSV. Keep lat/lon and
# the WorldCover class for later diagnostics and for the Phase 4 labeling
# sample.
#
# Cost note: these are batch table exports over global regions; run them in a
# month with EECU headroom.

library(AlphaSDM)

# ---- Config ------------------------------------------------------------------
YEAR         <- 2023                      # one recent year for v1
N_TOTAL      <- 3e6                       # within the plan's 2-5M range
DRIVE_FOLDER <- "alphasdm_sae_sample_v1"
SAMPLE_SCALE <- 10                        # native embedding resolution

# ESA WorldCover v200 classes. 80 (permanent water) is kept: the embeddings
# are masked over open water so those draws thin out on their own, and coastal
# water pixels that do carry embeddings are worth having.
WC_CLASSES <- c(10L, 20L, 30L, 40L, 50L, 60L, 70L, 80L, 90L, 95L, 100L)

AlphaSDM:::ensure_gee_authenticated()
ee <- reticulate::import("ee")

# Continents from LSIB simplified world regions, grouped into landmasses.
lsib <- ee$FeatureCollection("USDOS/LSIB_SIMPLE/2017")
continents <- list(
  africa        = c("AFRICA"),
  europe        = c("EUROPE"),
  asia          = c("E ASIA", "S ASIA", "SW ASIA", "CENTRAL ASIA"),
  north_america = c("NORTH AMERICA", "CENTRAL AMERICA", "CARIBBEAN"),
  south_america = c("SOUTH AMERICA"),
  oceania       = c("AUSTRALIA", "OCEANIA", "SE ASIA")
)

per_class_per_continent <- as.integer(ceiling(
  N_TOTAL / (length(WC_CLASSES) * length(continents))))
message(sprintf("Target: %d points per WorldCover class per continent (%d total).",
                per_class_per_continent,
                per_class_per_continent * length(WC_CLASSES) * length(continents)))

# The shared value convention: raw asset values. Everything the SAE ever sees
# must pass through alphaearth_rescale(), the same function the deployed
# encoder applies, so training and deployment can never diverge.
emb <- AlphaSDM:::alphaearth_rescale(AlphaSDM:::get_embedding_image(YEAR))
wc  <- ee$ImageCollection("ESA/WorldCover/v200")$first()$select("Map")$rename("wc_class")

stack <- emb$
  addBands(wc)$
  addBands(ee$Image$pixelLonLat())

for (cont in names(continents)) {
  region <- lsib$
    filter(ee$Filter$inList("wld_rgn", as.list(continents[[cont]])))$
    geometry()$
    dissolve(1000)

  sample <- stack$stratifiedSample(
    numPoints  = per_class_per_continent,
    classBand  = "wc_class",
    region     = region,
    scale      = SAMPLE_SCALE,
    classValues = as.list(WC_CLASSES),
    classPoints = as.list(rep(per_class_per_continent, length(WC_CLASSES))),
    seed       = 42L,
    geometries = FALSE,      # lon/lat travel as the pixelLonLat bands
    tileScale  = 16L
  )

  task <- ee$batch$Export$table$toDrive(
    collection     = sample,
    description    = sprintf("alphasdm_sae_sample_%s", cont),
    folder         = DRIVE_FOLDER,
    fileNamePrefix = sprintf("sae_sample_v1_%s", cont),
    fileFormat     = "CSV"
  )
  task$start()
  message(sprintf("Started export for %s.", cont))
}

# Record the sampling metadata beside the outputs. "rescale: raw" documents
# the value convention the sample was exported in.
meta <- list(year = YEAR, n_target = N_TOTAL, scale = SAMPLE_SCALE,
             stratification = "WorldCover v200 class x continent (LSIB)",
             rescale = "raw (alphaearth_rescale identity)",
             seed = 42L, classes = WC_CLASSES,
             started = format(Sys.time(), tz = "UTC"))
dir.create("dev/artifacts", showWarnings = FALSE, recursive = TRUE)
jsonlite::write_json(meta, "dev/artifacts/sae_sample_v1_metadata.json",
                     auto_unbox = TRUE, pretty = TRUE)
message("Monitor with sdm_gee_status(); download the Drive folder when done.")
