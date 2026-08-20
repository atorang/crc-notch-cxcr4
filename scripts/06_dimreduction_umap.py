"""06_dimreduction_umap.py

Inputs:  data/processed/pelka/Pelka_CXCR4_{high,low}_processed.h5ad
Outputs: results/figures/UMAP_CXCR4_high_low_patients.pdf
         results/figures/UMAP_3Doners.pdf
         data/processed/pelka/Pelka_CXCR4+_processed_{donor}.h5ad
"""
import sys
import time
from pathlib import Path

import numpy as np
import scanpy as sc
import yaml
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))
from crcnotch.dimred import run_dimred  # noqa: E402
from crcnotch.signature_scoring import load_gene_sets, score_signature_zscore  # noqa: E402
from crcnotch.umap_plots import plot_gene_or_score_panel, plot_patient_panel  # noqa: E402

cfg = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())
dimred_all = cfg["dimensionality_reduction"]
seed = cfg["seed"]

pelka_dir = REPO_ROOT / cfg["paths"]["data_processed"] / "pelka"
fig_dir = REPO_ROOT / cfg["paths"]["results_figures"]
fig_dir.mkdir(parents=True, exist_ok=True)
gene_sets = load_gene_sets(cfg, REPO_ROOT)

# ---- (a) combined CXCR4-high/low UMAP, colored by patient ------------------
adata_low = sc.read_h5ad(pelka_dir / "Pelka_CXCR4_low_processed.h5ad")
adata_high = sc.read_h5ad(pelka_dir / "Pelka_CXCR4_high_processed.h5ad")

adata_low.obsm["X_umap"][:, 0] += -35
adata_low.obsm["X_umap"][:, 1] += -5

adata_combined = sc.concat([adata_low, adata_high], join="inner")
adata_combined.X = adata_combined.raw.X.copy()
sc.pp.scale(adata_combined, max_value=cfg["dimensionality_reduction"]["cohort_level"]["scale_clip_max"])
adata_combined.obsm["X_umap"] = np.vstack([adata_low.obsm["X_umap"], adata_high.obsm["X_umap"]])
umap_xy = adata_combined.obsm["X_umap"]

fig, axes = plt.subplots(1, 2, figsize=(18, 4), constrained_layout=True)
plot_patient_panel(axes[0], fig, umap_xy, adata_combined.obs["donor_id"].values)
plot_gene_or_score_panel(axes[1], fig, adata_combined, umap_xy, "CXCR4", True, "#2D5A43")
out2 = fig_dir / "UMAP_CXCR4_high_low_patients.pdf"
fig.savefig(out2, bbox_inches="tight", dpi=300)
plt.close(fig)
print("wrote", out2)

# ---- (b) per-patient dimred for the 3-donor figure -------------------------
three_donors = dimred_all["three_donor_figure"]
per_patient_cfg = dimred_all["per_patient"]

adata_high_raw = adata_high  # already CXCR4-high, pre-scale (raw.X preserved)
for donor in three_donors:
    t0 = time.time()
    sub = adata_high_raw[adata_high_raw.obs["donor_id"] == donor].copy()
    sub.X = sub.raw.X.copy()  # reset to log-normalized before this donor's own HVG/scale/PCA
    sub = run_dimred(
        sub, n_hvg=per_patient_cfg["n_hvg"], n_pcs=per_patient_cfg["n_pcs_computed"],
        n_pcs_for_neighbors=per_patient_cfg["n_pcs_used"],
        n_neighbors=per_patient_cfg["n_neighbors"], leiden_resolution=per_patient_cfg["leiden_resolution"],
        seed=seed, scale_clip_max=dimred_all["cohort_level"]["scale_clip_max"],
        umap_min_dist=dimred_all["umap_min_dist"],
    )
    score_signature_zscore(sub, gene_sets, use_raw=True)
    out_path = pelka_dir / f"Pelka_CXCR4+_processed_{donor}.h5ad"
    sub.write_h5ad(out_path, compression="gzip")
    print(f"{donor}: {sub.shape} written to {out_path} ({time.time()-t0:.1f}s)")

# ---- 3-donor grid plot ------------------------------------------------------
row_panels = ["RSC", "CBC", "CXCR4", "Fetal", "wnt_msigdb"]
ncol, nrow = len(three_donors), len(row_panels)
fig, axes = plt.subplots(nrow, ncol, figsize=(4.5 * ncol, 3.5 * nrow), constrained_layout=True)
for r, panel in enumerate(row_panels):
    for c, donor in enumerate(three_donors):
        sub = sc.read_h5ad(pelka_dir / f"Pelka_CXCR4+_processed_{donor}.h5ad")
        ax = axes[r, c]
        base_color = "#7B328C" if panel in ("Fetal", "RSC") else "#2D5A43"
        plot_gene_or_score_panel(ax, fig, sub, sub.obsm["X_umap"], panel, panel == "CXCR4", base_color,
                                  point_size=25, donor=donor, title_fontsize=12)
out3 = fig_dir / "UMAP_3Doners.pdf"
fig.savefig(out3, bbox_inches="tight", dpi=300)
plt.close(fig)
print("wrote", out3)
