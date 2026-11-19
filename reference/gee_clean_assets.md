# Remove leftover AlphaSDM temporary assets

AlphaSDM writes temporary Earth Engine assets while it works and deletes
them when it finishes. A run that is killed or loses its connection
never reaches that cleanup, so the asset stays and counts against the
project's storage quota. This removes those leftovers.

## Usage

``` r
gee_clean_assets(older_than_hours = 48, dry_run = FALSE, project = NULL)
```

## Arguments

- older_than_hours:

  Keep assets younger than this. Default 48.

- dry_run:

  If TRUE, report what would be deleted and delete nothing.

- project:

  Earth Engine project id, or NULL to use the saved one.

## Value

The asset ids removed, or the ones that would be, invisibly.

## Details

Only assets named by AlphaSDM (\`alphasdm\_...\`) are considered. An
asset is kept if a task for it is still pending or running, or if it is
newer than \`older_than_hours\`, so a job in progress in another session
is not disturbed.

## Examples

``` r
if (FALSE) { # \dontrun{
gee_clean_assets(dry_run = TRUE)   # list what would be removed
} # }
```
