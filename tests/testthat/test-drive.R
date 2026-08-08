# The tile-name filter decides what gets deleted from the user's Drive, so its
# behaviour is pinned here. It must accept exactly what Earth Engine writes for an
# export and nothing a user might have named themselves.

test_that("names Earth Engine writes for an export are accepted", {
  p <- "alphasdm_20260807012507_ab12cd"
  expect_true(is_export_tile_name(p, paste0(p, "-0000016384-0000110592.tif")))
  expect_true(is_export_tile_name(p, paste0(p, "-0000000000-0000000000.tif")))
  # A region that fits in one file gets no tile suffix.
  expect_true(is_export_tile_name(p, paste0(p, ".tif")))
  # Multi-band exports insert a _rNNN band suffix before the tile suffix.
  expect_true(is_export_tile_name(p, paste0(p, "_r001-0000000000-0000004096.tif")))
  expect_true(is_export_tile_name(p, paste0(p, "_r017.tif")))
  # Grid-cell exports add a column suffix.
  expect_true(is_export_tile_name(p, paste0(p, "_r003_c012-0000000000-0000004096.tif")))
  expect_true(is_export_tile_name(p, paste0(p, "_r003_c012.tif")))
})

test_that("user files that merely contain the prefix are rejected", {
  p <- "alphasdm_20260807012507_ab12cd"
  expect_false(is_export_tile_name(p, paste0(p, "_notes.txt")))
  expect_false(is_export_tile_name(p, paste0(p, "-extra.tif")))
  expect_false(is_export_tile_name(p, paste0("backup_", p, ".tif")))
  expect_false(is_export_tile_name(p, "AlphaSDM QR code"))
  expect_false(is_export_tile_name(p, "AlphaSDM output"))
  # Same name, wrong extension.
  expect_false(is_export_tile_name(p, paste0(p, "-0000016384-0000110592.tiff")))
  # A band marker alone is not license for arbitrary trailing text.
  expect_false(is_export_tile_name(p, paste0(p, "_r001_backup.tif")))
  expect_false(is_export_tile_name(p, paste0(p, "_c012.tif")))
})

test_that("the filter is vectorised and keeps order", {
  p <- "x"
  nm <- c("x.tif", "y.tif", "x-1-2.tif", "x.tif.bak")
  expect_identical(is_export_tile_name(p, nm), c(TRUE, FALSE, TRUE, FALSE))
})
