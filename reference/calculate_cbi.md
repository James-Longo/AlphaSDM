# Calculate the Continuous Boyce Index

Measures how far the ratio of predicted to expected presences rises with
suitability. A well calibrated model gives a value near 1, a random one
near 0.

## Usage

``` r
calculate_cbi(pos_scores, all_scores, window_width = 0.1, n_bins = 100)
```

## Arguments

- pos_scores:

  Numeric suitability scores at presence points. NA is dropped.

- all_scores:

  Numeric suitability scores at all points, presence and background
  together. NA is dropped.

- window_width:

  Width of the moving window, as a proportion of the score range.

- n_bins:

  Number of window positions to evaluate.

## Value

A single number in \[-1, 1\], the Spearman correlation between window
position and the predicted-to-expected ratio. Returns 0 when there are
no presence scores, when all scores are equal, or when the correlation
is undefined.

## Examples

``` r
set.seed(1)
background <- runif(500)
presences  <- rbeta(50, 4, 2)   # presences sit at higher scores
calculate_cbi(presences, c(presences, background))
#> [1] 0.5549083
```
