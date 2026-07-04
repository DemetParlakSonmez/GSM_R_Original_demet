#!/usr/bin/env bash
# run-live.sh — keep the gsmR Shiny app running with live reload.
#
#   • Edit inst/shiny/app.R (UI/server)  → Shiny auto-reloads the browser
#     (shiny.autoreload=TRUE), no restart, no manual refresh.
#   • Edit R/*.R / DESCRIPTION / NAMESPACE (package code) → the R process is
#     restarted with a fresh devtools::load_all(); then just refresh the tab.
#
# Usage:
#   ./run-live.sh                 # runs on port 8760
#   ./run-live.sh 8899            # runs on a custom port
#   nohup ./run-live.sh >app-live.log 2>&1 &   # keep running after you log out
#
set -u
export PATH=/opt/conda/envs/gsmr/bin:$PATH
# Single-threaded BLAS so forked workers (BRIDGE mb/huge, glasso, MCCV) can't
# deadlock on OpenBLAS + fork(). Without this, BRIDGE "runs forever".
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 GOTO_NUM_THREADS=1
cd "$(dirname "$0")"
PORT="${1:-8760}"
APP_PID=""

free_port() {                     # wait until the port is actually released
  for _ in 1 2 3 4 5 6 7 8; do
    curl -s -o /dev/null "http://127.0.0.1:${PORT}/" && sleep 1 || return 0
  done
}

start() {
  free_port
  Rscript -e "options(shiny.autoreload=TRUE, gsmR.launch_dir=normalizePath('.'));
              suppressMessages(devtools::load_all('.'));
              suppressMessages({library(shiny);library(bslib);library(DT)});
              shiny::runApp('inst/shiny/app.R', port=${PORT}, host='0.0.0.0', launch.browser=FALSE)" &
  APP_PID=$!
  echo ">> app started (pid ${APP_PID}) — http://localhost:${PORT}"
}

stop() { [ -n "${APP_PID}" ] && { kill "${APP_PID}" 2>/dev/null; wait "${APP_PID}" 2>/dev/null; }; }
trap 'echo; echo ">> stopping"; stop; exit 0' INT TERM

start
echo ">> watching R/ DESCRIPTION NAMESPACE for changes (Ctrl-C to stop)"
while inotifywait -qq -e modify,create,delete,move -r R DESCRIPTION NAMESPACE; do
  echo ">> package code changed — restarting…"
  stop
  start
done
