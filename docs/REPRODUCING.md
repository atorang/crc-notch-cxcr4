# REPRODUCING.md

## Setup

```bash
conda env create -f env/environment.yml
conda activate crc-notch-cxcr4
```

## Running everything

```bash
bash run_all.sh          # Linux/macOS, or Git Bash on Windows
# or
.\run_all.ps1            # native PowerShell on Windows
# or, for partial reruns with dependency tracking:
snakemake -j1 --snakefile workflow/Snakefile all
```

## Determinism

Every script sets `seed: 42` (from `config.yaml`) for every stochastic step
(NumPy/Python `random`, PCA, neighbor graphs, UMAP, Leiden, Palantir,
K-means). 


