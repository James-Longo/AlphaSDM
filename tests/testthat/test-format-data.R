raw <- function(...) {
  d <- data.frame(
    lon  = c(-71.5, -71.6, -71.7),
    lat  = c(44.5, 44.6, 44.7),
    yr   = c(2020, 2021, 2022),
    occ  = c(1, 0, 1),
    stringsAsFactors = FALSE
  )
  modifyList(d, list(...))
}

test_that("columns are renamed to the package's standard names", {
  out <- format_data(raw(), coords = c("lon", "lat"), year = "yr", presence = "occ")
  expect_named(out, c("longitude", "latitude", "year", "present"))
  expect_equal(nrow(out), 3)
})

test_that("missing columns raise informative errors", {
  expect_error(format_data(raw(), coords = c("nope", "lat"), year = "yr"),
               "Coordinate columns not found")
  expect_error(format_data(raw(), coords = c("lon", "lat"), year = "nope"),
               "Year column not found")
  expect_error(format_data(raw(), coords = c("lon", "lat"), year = "yr", presence = "nope"),
               "Presence column not found")
  expect_error(format_data(raw(), coords = c("lon", "lat")), "argument \"year\" is missing")
})

test_that("records outside Alpha Earth coverage are dropped", {
  d <- raw(yr = c(2005, 2021, 2022))
  expect_message(out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ"),
                 "temporal coverage")
  expect_equal(nrow(out), 2)
  expect_true(all(out$year >= 2017 & out$year <= 2025))
})

test_that("date columns are reduced to their year", {
  # All presences, so the presence-first reordering below cannot shuffle the rows.
  d <- raw(yr = as.Date(c("2020-05-01", "2021-06-15", "2022-07-30")), occ = c(1, 1, 1))
  out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ")
  expect_equal(out$year, c(2020, 2021, 2022))
})

test_that("omitting presence marks every record as a presence", {
  expect_message(out <- format_data(raw(), coords = c("lon", "lat"), year = "yr"),
                 "presence-only")
  expect_true(all(out$present == 1))
})

test_that("duplicate coordinates collapse presence-first", {
  # Same coordinate and year recorded as both an absence and a presence: the
  # presence must win, and it must be the row that survives.
  d <- data.frame(lon = c(-71.5, -71.5), lat = c(44.5, 44.5),
                  yr = c(2020, 2020), occ = c(0, 1))
  expect_message(out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ"),
                 "duplicate or conflicting")
  expect_equal(nrow(out), 1)
  expect_equal(out$present, 1)
})

test_that("presences are ordered first (the libsvm row-order contract)", {
  # GEE's libsvm takes its positive class from the first label it sees, and the
  # SVM score flip downstream assumes that class is presence.
  d <- data.frame(lon = c(-71.5, -71.6, -71.7), lat = c(44.5, 44.6, 44.7),
                  yr = rep(2020, 3), occ = c(0, 1, 0))
  out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ")
  expect_equal(out$present[1], 1)
})

test_that("rows with missing values are removed", {
  d <- raw(lat = c(44.5, NA, 44.7))
  expect_message(out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ"),
                 "missing \\(NA\\) values")
  expect_equal(nrow(out), 2)
})

test_that("swapped coordinates are flagged", {
  d <- data.frame(lon = c(44.5, 44.6), lat = c(-171.5, -171.6),
                  yr = c(2020, 2021), occ = c(1, 1))
  expect_warning(format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ"),
                 "SWAPPED")
})

test_that("a species column is carried through when requested", {
  d <- cbind(raw(), sp = c("a", "b", "a"))
  out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ", species = "sp")
  expect_true("species" %in% names(out))
})
