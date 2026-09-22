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
