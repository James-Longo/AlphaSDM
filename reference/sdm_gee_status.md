# Show the status of recent AlphaSDM Earth Engine batch tasks

Lists the export and classifier tasks AlphaSDM has started, with each
task's state and age. Call it from a second session to see what Earth
Engine is doing while a long run is in progress.

## Usage

``` r
sdm_gee_status(active_only = TRUE, since_minutes = 180)
```

## Arguments

- active_only:

  TRUE shows only pending and running tasks. FALSE also lists recently
  finished ones.

- since_minutes:

  Only include tasks created within this many minutes.

## Value

A data frame of tasks with \`description\`, \`state\` and \`age_min\`,
invisibly. Also prints them.
