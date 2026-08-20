# crc-notch-cxcr4

Reproducible analysis code for the manuscript *"CXCR4 dependent dynamic switching between WNT- and YAP-active cell populations drives aggressive behavior in colorectal cancer"* 

## Summary

This repository contains the computational pipeline and statistical code supporting the analysis of single-cell transcriptomics, and bulk RNA-seq data of colorectal cancer.


## Quick start

```bash
conda env create -f env/environment.yml && conda activate crc-notch-cxcr4
bash run_all.sh            # or .\run_all.ps1 on Windows, or `snakemake -j1 --snakefile workflow/Snakefile all`
```

See `docs/REPRODUCING.md` for per-step runtime/memory and determinism
notes, and `docs/DATA_ACCESS.md` for what's automated vs. requires a manual
download step (Pelka/SCP1162).

## Repository layout

`config/config.yaml`: All parameters, thresholds, and random seeds used in the pipeline.

`R/`, `python/crcnotch/`: Core utility functions and standard implementations for signature scoring, GSVA, trajectory and survival analysis.

`scripts/00`–`10`: Numbered, sequential pipeline execution scripts.

## Citation
