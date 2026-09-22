# Calculate discrimination and calibration metrics

Calculate discrimination and calibration metrics

## Usage

``` r
calculate_classifier_metrics(scores_pos, scores_neg)
```

## Arguments

- scores_pos:

  Numeric suitability scores at presence points. NA is dropped.

- scores_neg:

  Numeric suitability scores at background or absence points. NA is
  dropped.

## Value

A named list: \`cbi\`, \`auc_roc\`, \`auc_prg\`, \`tss\`, \`ba\` and
\`cor\`, each a single number. Scores need only be on a common scale
within one call, since every metric except \`cor\` depends on the
ranking alone. When either class is empty the list is filled with the
no-skill values.

## Examples

``` r
set.seed(1)
presences <- rbeta(50, 4, 2)
absences  <- rbeta(200, 2, 4)
str(calculate_classifier_metrics(presences, absences))
#> List of 6
#>  $ cbi    : num 0.988
#>  $ auc_roc: num 0.895
#>  $ auc_prg: num 0.908
#>  $ tss    : num 0.64
#>  $ ba     : num 0.82
#>  $ cor    : num 0.59
```
