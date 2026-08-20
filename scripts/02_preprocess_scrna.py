"""02_preprocess_scrna.py

QC-filters the Pelka et al. malignant-epithelial scRNA-seq compartment:
extracts the paper's `pEpi*` metagene columns into `.obs`, then applies the
QC thresholds from config.yaml (min/max genes, max %MT, min cells/gene).

Input:  data/raw/scrnaseq/scRNAseq_Pelka_Epithelial.h5ad
Output: data/processed/pelka/pelka_epithelial_qc.h5ad
"""
import sys
import time
from pathlib import Path

import scanpy as sc
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))
from crcnotch.qc import extract_pelka_metaprograms, qc_filter  # noqa: E402

cfg = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())

raw_path = REPO_ROOT / cfg["paths"]["data_raw"] / "scrnaseq" / "scRNAseq_Pelka_Epithelial.h5ad"
out_dir = REPO_ROOT / cfg["paths"]["data_processed"] / "pelka"
out_dir.mkdir(parents=True, exist_ok=True)
out_path = out_dir / "pelka_epithelial_qc.h5ad"

t0 = time.time()
adata = sc.read_h5ad(raw_path)
adata = adata[adata.obs["ClusterMidway"] == "EpiT"].copy()
adata.var["ensemble_id"] = adata.var_names
adata.var_names = adata.var["gene_symbol"].values.astype(str)
adata.var_names_make_unique()

adata = extract_pelka_metaprograms(adata)

qc = cfg["qc_scrna"]
adata = qc_filter(
    adata,
    min_genes=qc["min_genes"], max_genes=qc["max_genes"],
    max_pct_mt=qc["max_pct_mt"], min_cells_per_gene=qc["min_cells_per_gene"],
)

print(f"QC-filtered: {adata.shape[0]} cells x {adata.shape[1]} genes ({time.time()-t0:.1f}s)")
adata.write_h5ad(out_path, compression="gzip")
print("wrote", out_path)
