# Remove leftover AlphaSDM temporary assets

The package writes temporary Earth Engine assets while it works and
deletes them when it finishes. A run that is killed, crashes, or loses
its connection never reaches that cleanup, so the asset stays and counts
against the project's storage quota. This removes those leftovers.

## Usage

``` r
sdm_clean_assets(
  older_than_hours = 48,
  dry_run = FALSE,
  project = NULL,
  quiet = FALSE
)
```

## Arguments

- older_than_hours:

  Keep assets younger than this. Default 48.

- dry_run:

  If TRUE, report what would be deleted and delete nothing.

- project:

  Earth Engine project id, or NULL to use the saved one.

- quiet:

  If TRUE, do not print anything.

## Value

The asset ids removed, or the ones that would be, invisibly.

## Details

Only assets this package created are considered, matched on the
\`alphasdm\_\` name it gives them. An asset is kept if a batch task for
it is still pending or running, or if it is newer than
\`older_than_hours\`, so a job in progress in another session is not
disturbed. A large export can run for many hours, which is why the
default is deliberately generous.
