Write-Host "Starting GSM-R Shiny Application..." -ForegroundColor Cyan

# Define the path to Rscript.exe
$RscriptPath = "C:\Program Files\R\R-4.4.3\bin\Rscript.exe"

# Check if R is installed at the specified path
if (-Not (Test-Path $RscriptPath)) {
    Write-Host "ERROR: Rscript.exe not found! Please ensure R-4.4.3 is installed in C:\Program Files\R\" -ForegroundColor Red
    Write-Host "If you are using a different version, open start-app.ps1 and update the RscriptPath variable." -ForegroundColor Yellow
    exit
}

Write-Host "Checking required packages and opening the application in browser at 127.0.0.1:8760..." -ForegroundColor Green

# Launch the R application
& $RscriptPath -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); pkgs <- c('pkgload', 'shiny', 'bslib', 'DT', 'data.table', 'randomForest', 'e1071', 'huge'); new_pkgs <- pkgs[!(pkgs %in% installed.packages()[,'Package'])]; if(length(new_pkgs)) install.packages(new_pkgs); pkgload::load_all('.'); shiny::runApp('inst/shiny/app.R', port = 8760, host = '127.0.0.1', launch.browser = TRUE)"
