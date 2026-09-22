# Turn AlphaSDM console output on or off

Every progress line the package prints is a message, not a warning or a
print, so this suppresses all of them at once. Errors are unaffected.

## Usage

``` r
sdm_verbose(verbose = TRUE)
```

## Arguments

- verbose:

  Logical. FALSE suppresses all non-error output.

## Value

\`verbose\`, invisibly. Called for its effect on package state.

## Examples

``` r
sdm_verbose(FALSE)   # silence progress messages
sdm_verbose(TRUE)
```
