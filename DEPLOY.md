# Deploying the gsmR dashboard (always-on)

The gsmR Shiny GUI normally runs locally via `gsmR::run_app()`, which only
lives as long as your R session. To host it **always-on** with a permanent
URL, use one of the options below.

## Option A — Docker (recommended, portable)

A [`Dockerfile`](Dockerfile) is included that bakes in R, every package
dependency, and the gsmR package itself, then serves the dashboard on port
8080.

```bash
docker build -t gsmr .
docker run -p 8080:8080 gsmr
# open http://localhost:8080
```

This same image runs unchanged on any container host:

| Host | Notes |
|------|-------|
| **Hugging Face Spaces** | Create a Space with the *Docker* SDK, push this repo. Free CPU — good fit for the parallel MCCV workload. |
| **Google Cloud Run** | `gcloud run deploy --source .` — scales to zero, pay per request. |
| **Render / Railway / Fly.io** | Point at this repo; they detect the `Dockerfile` automatically. |
| **Self-hosted VM** | `docker run` on any DigitalOcean/EC2 box for 24/7 hosting with no usage caps. |

## Option B — shinyapps.io (no Docker, quickest URL)

Posit's hosted Shiny service. Free tier: 5 apps, 25 active hours/month.

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "<acct>", token = "<token>", secret = "<secret>")
rsconnect::deployApp(appDir = "inst/shiny", appName = "gsmR")
```

Because the app depends on the `gsmR` package being available on the server,
publish this repo to GitHub and let the app install it on startup, or vendor
the package functions alongside `app.R`.

## Local development

```r
library(gsmR)
run_app()   # opens in your browser, stops when the session ends
```
