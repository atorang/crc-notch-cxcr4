"""04_signature_scores_scrna.py

Input:  data/processed/pelka/Pelka_CXCR4_high_processed.h5ad
Output: results/tables/signature_scores_scrna_cxcr4_high.tsv
"""
import sys
from pathlib import Path

import scanpy as sc
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))
from crcnotch.signature_scoring import load_gene_sets, score_signature_zscore  # noqa: E402

cfg = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())

pelka_dir = REPO_ROOT / cfg["paths"]["data_processed"] / "pelka"
table_dir = REPO_ROOT / cfg["paths"]["results_tables"]
table_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(pelka_dir / "Pelka_CXCR4_high_processed.h5ad")
gene_sets = load_gene_sets(cfg, REPO_ROOT)
score_signature_zscore(adata, gene_sets, use_raw=True)

score_cols = list(gene_sets.keys())
out_path = table_dir / "signature_scores_scrna_cxcr4_high.tsv"
adata.obs[["donor_id"] + score_cols].to_csv(out_path, sep="\t")
print("wrote", out_path)
