# Forget the saved Earth Engine project, and optionally sign out

Removes the project ID AlphaSDM saved. In an interactive session it then
offers to delete the Earth Engine sign-in credentials as well; those are
shared by every tool on this computer that uses Earth Engine, so they
are kept unless you agree. Run
[`setup_gee`](https://james-longo.github.io/AlphaSDM/reference/setup_gee.md)
afterwards to reconnect.

## Usage

``` r
clear_gee_credentials()
```

## Value

Invisibly, `TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
clear_gee_credentials()
} # }
```
