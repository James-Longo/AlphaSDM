# Show recent AlphaSDM Earth Engine tasks

Lists the export tasks AlphaSDM has started, with each task's state and
age. Call it from a second R session to see what Earth Engine is doing
while a long run is in progress. To check the connection itself, use
\[gee_status()\].

## Usage

``` r
gee_tasks(active_only = TRUE, since_minutes = 180)
```

## Arguments

- active_only:

  TRUE shows only pending and running tasks; FALSE also lists recently
  finished ones.

- since_minutes:

  Only include tasks created within this many minutes.

## Value

A data frame of tasks with \`description\`, \`state\` and \`age_min\`,
invisibly. Also prints them.

## Examples

``` r
if (FALSE) { # \dontrun{
gee_tasks(active_only = FALSE)
} # }
```
