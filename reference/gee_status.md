# Report the Google Earth Engine connection status

Prints whether the Earth Engine client is available, whether sign-in
credentials exist and of which kind, which project is configured, and
whether a live connection succeeds. To monitor running Earth Engine
tasks, use
[`sdm_gee_status`](https://james-longo.github.io/AlphaSDM/reference/sdm_gee_status.md)
instead.

## Usage

``` r
gee_status(check_live = TRUE)
```

## Arguments

- check_live:

  If `TRUE` (default), make a small request to confirm that the
  credentials work, not only that they are on disk.

## Value

Invisibly, a named list of the status fields.

## Examples

``` r
if (FALSE) { # \dontrun{
gee_status()
} # }
```
