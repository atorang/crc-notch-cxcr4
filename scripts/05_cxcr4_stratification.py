"""05_cxcr4_stratification.py

Input:  data/processed/pelka/pelka_epithelial_qc.h5ad
Output: results/figures/Barplots_CXCR4_pos_neg.png
        data/processed/pelka/{Pelka_CXCR4_high,Pelka_CXCR4_low}_processed.h5ad
"""
import sys
import time
from pathlib import Path

import scanpy as sc
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))
from crcnotch.cxcr4_stratification import (  # noqa: E402
    compute_cxcr4_positivity, per_patient_fraction, plot_ranked_fraction_barplot, stratify_patients,
)
from crcnotch.dimred import run_dimred  # noqa: E402
from crcnotch.qc import ensure_log_normalized  # noqa: E402
from crcnotch.signature_scoring import load_gene_sets, score_signature_zscore  # noqa: E402

cfg = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())
qc, strat, dimred_cfg = cfg["qc_scrna"], cfg["cxcr4_stratification"], cfg["dimensionality_reduction"]["cohort_level"]
seed = cfg["seed"]

in_path = REPO_ROOT / cfg["paths"]["data_processed"] / "pelka" / "pelka_epithelial_qc.h5ad"
out_dir = REPO_ROOT / cfg["paths"]["data_processed"] / "pelka"
fig_dir = REPO_ROOT / cfg["paths"]["results_figures"]
fig_dir.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(in_path)
adata = ensure_log_normalized(adata)
adata = compute_cxcr4_positivity(adata, positivity_threshold=strat["positivity_threshold"])

per_pat = per_patient_fraction(adata, patient_key="donor_id", min_tumor_cells=qc["min_cells_per_tumor"])
per_pat = stratify_patients(per_pat, strat["high_fraction_cutoff"], strat["low_fraction_cutoff"])

barplot_path = fig_dir / "Barplots_CXCR4_pos_neg.png"
plot_ranked_fraction_barplot(per_pat, strat["high_fraction_cutoff"], strat["low_fraction_cutoff"], barplot_path)
print("wrote", barplot_path)
print(per_pat["CXCR4_group"].value_counts())

grp_map = per_pat["CXCR4_group"].to_dict()
adata.obs["CXCR4_group"] = adata.obs["donor_id"].map(grp_map)
adata = adata[adata.obs["CXCR4_group"].notna()].copy()
print("tumor cells used in workflow:", adata.shape)
print(adata.obs["CXCR4_group"].value_counts())

gene_sets = load_gene_sets(cfg, REPO_ROOT)

for group_label, out_name in [("CXCR4_high", "Pelka_CXCR4_high_processed.h5ad"),
                                ("CXCR4_low", "Pelka_CXCR4_low_processed.h5ad")]:
    t0 = time.time()
    sub = adata[adata.obs["CXCR4_group"] == group_label].copy()
    sub = run_dimred(
        sub, n_hvg=dimred_cfg["n_hvg"], n_pcs=dimred_cfg["n_pcs"], n_neighbors=dimred_cfg["n_neighbors"],
        leiden_resolution=dimred_cfg["leiden_resolution"], seed=seed,
        hvg_flavor=dimred_cfg["hvg_flavor"], scale_clip_max=dimred_cfg["scale_clip_max"],
        leiden_flavor=dimred_cfg["leiden_flavor"], leiden_n_iterations=dimred_cfg["leiden_n_iterations"],
        pca_solver=dimred_cfg["pca_solver"],
    )
    score_signature_zscore(sub, gene_sets, use_raw=True)
    out_path = out_dir / out_name
    sub.write_h5ad(out_path, compression="gzip")
    print(f"{group_label}: {sub.shape} written to {out_path} ({time.time()-t0:.1f}s)")
