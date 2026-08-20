# Runs the full pipeline in order, stopping on the first error.
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

$RSCRIPT = if ($env:RSCRIPT) { $env:RSCRIPT } else { "Rscript" }
$PYTHON = if ($env:PYTHON) { $env:PYTHON } else { "python" }

function Run($cmd, $args_) {
    & $cmd @args_
    if ($LASTEXITCODE -ne 0) { throw "$cmd $($args_ -join ' ') failed with exit code $LASTEXITCODE" }
}

Write-Host "== 00: download data =="
Run $RSCRIPT @("scripts/00_download_data.R")
Run $PYTHON @("scripts/00_download_data.py")

Write-Host "== 01: preprocess bulk =="
Run $RSCRIPT @("scripts/01_preprocess_bulk.R")

Write-Host "== 02: preprocess scRNA =="
Run $PYTHON @("scripts/02_preprocess_scrna.py")

Write-Host "== 03: signature scores (bulk) =="
Run $RSCRIPT @("scripts/03_signature_scores_bulk.R")

Write-Host "== 04: signature scores (scRNA) =="
Run $PYTHON @("scripts/04_signature_scores_scrna.py")

Write-Host "== 05: CXCR4 stratification =="
Run $PYTHON @("scripts/05_cxcr4_stratification.py")

Write-Host "== 06: dimensionality reduction / UMAP =="
Run $PYTHON @("scripts/06_dimreduction_umap.py")

Write-Host "== 07: Palantir trajectory =="
Run $PYTHON @("scripts/07_trajectory_palantir.py")

Write-Host "== 08: survival (Kaplan-Meier) =="
Run $RSCRIPT @("scripts/08_survival_km.R")

Write-Host "== 09: mutation/expression boxplots =="
Run $RSCRIPT @("scripts/09_mutation_expression_boxplots.R")

Write-Host "== 10: assemble figures =="
Run $PYTHON @("scripts/10_assemble_figures.py")

Write-Host "Done. See results/figures/, results/tables/, results/logs/."
