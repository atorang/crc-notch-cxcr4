#!/usr/bin/env bash
# Runs the full pipeline in order, stopping on the first error.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

RSCRIPT="${RSCRIPT:-Rscript}"
PYTHON="${PYTHON:-python}"

echo "== 00: download data =="
"$RSCRIPT" scripts/00_download_data.R
"$PYTHON" scripts/00_download_data.py

echo "== 01: preprocess bulk =="
"$RSCRIPT" scripts/01_preprocess_bulk.R

echo "== 02: preprocess scRNA =="
"$PYTHON" scripts/02_preprocess_scrna.py

echo "== 03: signature scores (bulk) =="
"$RSCRIPT" scripts/03_signature_scores_bulk.R

echo "== 04: signature scores (scRNA) =="
"$PYTHON" scripts/04_signature_scores_scrna.py

echo "== 05: CXCR4 stratification =="
"$PYTHON" scripts/05_cxcr4_stratification.py

echo "== 06: dimensionality reduction / UMAP =="
"$PYTHON" scripts/06_dimreduction_umap.py

echo "== 07: Palantir trajectory =="
"$PYTHON" scripts/07_trajectory_palantir.py

echo "== 08: survival (Kaplan-Meier) =="
"$RSCRIPT" scripts/08_survival_km.R

echo "== 09: mutation/expression boxplots =="
"$RSCRIPT" scripts/09_mutation_expression_boxplots.R

echo "== 10: assemble figures =="
"$PYTHON" scripts/10_assemble_figures.py

echo "Done. See results/figures/, results/tables/, results/logs/."
