pkgload::load_all(".", quiet = TRUE)
options(shiny.port = 8760, shiny.host = "127.0.0.1")
cat("Starting Shiny app at http://127.0.0.1:8760 ...\n")
shiny::runApp("inst/shiny/app.R", port = 8760, host = "127.0.0.1", launch.browser = FALSE)
