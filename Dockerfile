# gsmR — always-on Shiny dashboard
#
# Builds a self-contained image with R, all package dependencies, and the
# gsmR package installed, then serves the Shiny GUI on port 8080.
#
#   docker build -t gsmr .
#   docker run -p 8080:8080 gsmr
#
# Then open http://localhost:8080
#
# This image is also directly deployable to Hugging Face Spaces (Docker SDK),
# Google Cloud Run, Render, Railway, Fly.io, or any host that runs containers.

FROM rocker/r-ver:4.3.3

# System libraries required by the R package stack (igraph, openxlsx, curl, xml).
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libglpk-dev \
        libgmp-dev \
        zlib1g-dev \
        libuv1 \
    && rm -rf /var/lib/apt/lists/*

# Use Posit Package Manager binaries for fast, reproducible installs.
RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))' \
        >> /usr/local/lib/R/etc/Rprofile.site

# Runtime + Suggests packages the Shiny app and full G-S-M workflow use.
RUN R -e 'install.packages(c( \
        "shiny", "bslib", "DT", "data.table", "stats", \
        "randomForest", "e1071", "caret", "igraph", "Matrix", \
        "RobustRankAggreg", "openxlsx", "remotes" \
    ))'

# Additional packages the Shiny GUI loads at startup / on demand:
#   visNetwork — interactive detection-network graph (referenced in the UI)
#   yaml       — settings import/export (a hard Import of gsmR)
#   xgboost    — optional XGBoost classifier
RUN R -e 'install.packages(c("visNetwork", "yaml", "xgboost"))'

WORKDIR /app
COPY . /app

# Install the gsmR package itself from the copied source.
RUN R -e 'remotes::install_local(".", dependencies = FALSE, upgrade = "never")'

EXPOSE 8080

# Bind to all interfaces so the container is reachable from the host.
CMD ["R", "-e", "options(gsmR.launch_dir = getwd()); library(gsmR); shiny::runApp(system.file('shiny', package = 'gsmR'), host = '0.0.0.0', port = 8080, launch.browser = FALSE)"]
