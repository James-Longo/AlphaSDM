# Report the Google Earth Engine Connection Status

Prints a quick diagnostic of how AlphaSDM is connected to Earth Engine:
which Python environment is bound, whether saved credentials exist and
are the personal-account (OAuth) type, which project is configured, and
whether a live connection succeeds. Use it to confirm setup or to
troubleshoot.

## Usage

``` r
gee_status(check_live = TRUE)
```

## Arguments

- check_live:

  If `TRUE` (default), perform a small server round-trip to confirm the
  credentials actually work, not just that they are on disk.

## Value

Invisibly, a named list of the status fields.

## Details

Note: this reports the \*connection\*. To monitor running server-side
export tasks, use
[`sdm_gee_status`](https://james-longo.github.io/AlphaSDM/reference/sdm_gee_status.md)
instead.
