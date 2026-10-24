# Rule: "run the app" (Uygulamayı Başlatma Kuralı)

When the user asks to "run the app", "uygulamayı başlat", or "run the app ile çalıştır", you must NOT search for the `Rscript.exe` path using `Get-ChildItem` (which causes hanging). Instead, use the `run_command` tool to execute the exact command below in the background (`IsDaemon: true`).

**Command to run:**
```powershell
& "C:\Program Files\R\R-4.4.3\bin\Rscript.exe" -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); pkgs <- c('pkgload', 'shiny', 'bslib', 'DT', 'data.table', 'randomForest', 'e1071'); new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,'Package'])]; if(length(new_pkgs)) install.packages(new_pkgs); pkgload::load_all('.'); shiny::runApp('inst/shiny/app.R', port = 8760, host = '127.0.0.1', launch.browser = TRUE)"
```

**Instructions for Agent:**
1. Do not use `Get-Command` or `Get-ChildItem` to locate Rscript. 
2. If the user's R version is different, prompt the user to specify their exact R version path rather than scanning the entire `C:\Program Files\R` directory recursively.
3. Always launch it as a background task (`IsDaemon: true`) and inform the user that it will be available at `http://127.0.0.1:8760`.
