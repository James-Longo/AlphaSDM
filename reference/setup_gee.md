# Set Up Google Earth Engine for AlphaSDM (one-time)

Connects AlphaSDM to Google Earth Engine using your personal Google
account. You only ever need to run this once per machine: it
authenticates through a single browser click and saves long-lived
credentials, so every future R session connects automatically with no
further prompts.

## Usage

``` r
setup_gee(project = NULL, force = FALSE, auth_mode = NULL)
```

## Arguments

- project:

  Google Cloud / Earth Engine project ID (e.g., `"my-ee-project"`). If
  `NULL`, the saved project is reused, or you are prompted
  interactively.

- force:

  If `TRUE`, re-authenticate even if valid credentials already exist.

- auth_mode:

  Optional override for the Earth Engine authorization flow, passed
  straight to `ee$Authenticate()`. Leave `NULL` (recommended) to let
  AlphaSDM choose: `"localhost"` (one-click, no paste) on a machine with
  a browser, and `"notebook"` on a detected headless/remote session. Set
  `"notebook"` yourself to force the paste-a-code flow, or `"gcloud"` if
  you use the gcloud CLI.

## Value

Invisibly `TRUE` on success, or `FALSE` if a step (such as the Python
install) requires you to restart R and re-run.

## Details

**Before you start** you need a free Earth Engine account. Sign up at
<https://earthengine.google.com/signup/>. Earth Engine is free for
noncommercial, research, education, and nonprofit use. Registration
links your Google account to a Cloud project (its ID is what you pass as
`project`).

**The authentication is a browser click, not a code to paste.** On a
desktop or laptop, `setup_gee()` opens your browser, you click
**Allow**, and the credential is captured automatically over a local
loopback port (`auth_mode = "localhost"`). Nothing is copied or pasted,
and the saved credentials do not expire with normal use.

Re-running `setup_gee()` when you are already connected is a harmless
no-op; it detects the working credentials and returns immediately.
