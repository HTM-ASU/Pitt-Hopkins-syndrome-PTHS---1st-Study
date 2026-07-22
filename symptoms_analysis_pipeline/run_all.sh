#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

mkdir -p output/tables output/figures output/logs data/processed

echo "[1/3] Preparing change scores..."
python scripts/python/01_prepare_change_scores.py | tee output/logs/01_prepare_change_scores.log

echo "[2/3] Running t-tests..."
python scripts/python/02_run_ttests.py | tee output/logs/02_run_ttests.log

echo "[3/3] Generating plots..."
Rscript scripts/R/03_plot_results.R | tee output/logs/03_plot_results.log

echo "Done. Results are in output/tables and output/figures."
