# Connect AlphaSDM to Google Earth Engine (one-time)

Signs in to Google Earth Engine with your own Google account and saves
your project ID, so later R sessions connect on their own. Run it once
per machine; running it again when already connected does nothing.

## Usage

``` r
setup_gee(project = NULL, force = FALSE, auth_mode = NULL)
```

## Arguments

- project:

  Earth Engine (Google Cloud) project ID, for example `"my-ee-project"`.
  If `NULL`, the saved project is used, or you are asked for one.

- force:

  If `TRUE`, sign in again even if valid credentials exist.

- auth_mode:

  Earth Engine sign-in flow, passed to `ee.Authenticate()`. Leave `NULL`
  to use `"localhost"` (a browser click) where a browser is available
  and `"notebook"` (paste a code) where it is not. `"gcloud"` uses the
  gcloud command-line tool.

## Value

Invisibly, `TRUE` once connected.

## Details

You need a free Earth Engine account first: register at
<https://earthengine.google.com/signup/>. Earth Engine is free for
noncommercial, research, education and nonprofit use, and registration
gives you the Cloud project ID to pass as `project`.

Signing in opens your browser, where you click **Allow**. The Earth
Engine client stores the resulting credentials in its own configuration
folder, as it does for every tool that uses Earth Engine. AlphaSDM saves
only the project ID, in `tools::R_user_dir("AlphaSDM", "config")`.

The Earth Engine Python client (`earthengine-api`) is provided through
reticulate, which sets up a Python environment for it the first time it
is needed, unless you have pointed reticulate at a Python of your own.

## Examples

``` r
if (FALSE) { # \dontrun{
# Needs an Earth Engine account and an interactive session.
setup_gee(project = "my-ee-project")
} # }
```
