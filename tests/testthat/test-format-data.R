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
  local_coverage_window()
  out <- format_data(raw(), coords = c("lon", "lat"), year = "yr", presence = "occ")
  expect_named(out, c("longitude", "latitude", "year", "present"))
  expect_equal(nrow(out), 3)
})

test_that("missing columns raise informative errors", {
  local_coverage_window()
  expect_error(format_data(raw(), coords = c("nope", "lat"), year = "yr"),
               "Coordinate columns not found")
  expect_error(format_data(raw(), coords = c("lon", "lat"), year = "nope"),
               "Year column not found")
  expect_error(format_data(raw(), coords = c("lon", "lat"), year = "yr", presence = "nope"),
               "Presence column not found")
  expect_error(format_data(raw(), coords = c("lon", "lat")), "argument \"year\" is missing")
})

test_that("records outside Alpha Earth coverage are dropped", {
  local_coverage_window()
  d <- raw(yr = c(2005, 2021, 2022))
  expect_message(out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ"),
                 "temporal coverage")
  expect_equal(nrow(out), 2)
  expect_true(all(out$year >= 2017 & out$year <= 2025))
})

test_that("date columns are reduced to their year", {
  local_coverage_window()
  # All presences, so format_data's presence-first sort cannot reorder the rows.
  d <- raw(yr = as.Date(c("2020-05-01", "2021-06-15", "2022-07-30")), occ = c(1, 1, 1))
  out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ")
  expect_equal(out$year, c(2020, 2021, 2022))
})

test_that("omitting presence marks every record as a presence", {
  local_coverage_window()
  expect_message(out <- format_data(raw(), coords = c("lon", "lat"), year = "yr"),
                 "presence-only")
  expect_true(all(out$present == 1))
})

test_that("duplicate coordinates collapse presence-first", {
  local_coverage_window()
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
  local_coverage_window()
  # GEE's libsvm takes its positive class from the first label it sees, and the
  # SVM score flip downstream assumes that class is presence.
  d <- data.frame(lon = c(-71.5, -71.6, -71.7), lat = c(44.5, 44.6, 44.7),
                  yr = rep(2020, 3), occ = c(0, 1, 0))
  out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ")
  expect_equal(out$present[1], 1)
})

test_that("rows with missing values are removed", {
  local_coverage_window()
  d <- raw(lat = c(44.5, NA, 44.7))
  expect_message(out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ"),
                 "missing \\(NA\\) values")
  expect_equal(nrow(out), 2)
})

test_that("swapped coordinates are rejected with a swap hint", {
  local_coverage_window()
  d <- data.frame(lon = c(44.5, 44.6), lat = c(-171.5, -171.6),
                  yr = c(2020, 2021), occ = c(1, 1))
  expect_error(format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ"),
               "SWAPPED")
})

test_that("a species column is carried through when requested", {
  local_coverage_window()
  d <- cbind(raw(), sp = c("a", "b", "a"))
  out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ", species = "sp")
  expect_true("species" %in% names(out))
})

test_that("an all-NA coordinate column is handled, not fatal", {
  local_coverage_window()
  # x[1] on an empty vector is NA, so the swapped-coordinate guard has to test for
  # NA rather than NULL or the comparison below it errors.
  d <- data.frame(lon = c(NA_real_, NA_real_), lat = c(NA_real_, NA_real_),
                  yr = c(2020, 2021), occ = c(1, 1))
  expect_message(out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ"),
                 "missing \\(NA\\) values")
  expect_equal(nrow(out), 0)
})

test_that("the coverage window comes from Earth Engine and is cached", {
  local_coverage_window(c(2019L, 2021L))
  expect_equal(alphaearth_year_range(), c(2019L, 2021L))   # cache hit, no call made
})

test_that("the reported window drives the filter and the advice", {
  local_coverage_window(c(2019L, 2021L))
  d <- data.frame(lon = rep(-71.5, 4), lat = rep(44.5, 4),
                  yr = c(2018, 2019, 2021, 2022), occ = rep(1, 4))
  expect_message(out <- format_data(d, coords = c("lon", "lat"), year = "yr", presence = "occ"),
                 "2019-2021")
  expect_true(all(out$year >= 2019 & out$year <= 2021))
})

test_that("projected coordinates are rejected with a WGS84 message", {
  local_coverage_window()
  utm <- data.frame(x = c(500000, 501000), y = c(4900000, 4901000),
                    yr = 2020L, occ = 1L)
  expect_error(
    format_data(utm, coords = c("x", "y"), year = "yr", presence = "occ"),
    "WGS84")
})

test_that("valid decimal degrees are not caught by the projection guard", {
  local_coverage_window()
  ok <- data.frame(lon = c(-68.5, 179.9, -179.9), lat = c(44.7, -89.9, 89.9),
                   yr = 2020L, occ = 1L)
  expect_silent(suppressMessages(
    format_data(ok, coords = c("lon", "lat"), year = "yr", presence = "occ")))
})
