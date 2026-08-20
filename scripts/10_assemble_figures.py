"""10_assemble_figures.py

"""
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())
fig_dir = REPO_ROOT / cfg["paths"]["results_figures"]

EXPECTED = {
    "KaplanMeier_TCGA_Nunes.pdf": "08_survival_km.R — TCGA+Nunes NOTCH-quartile survival",
    "HES1_HES5_by_FBXW7_WNT_Nunes.pdf": "09_mutation_expression_boxplots.R — Nunes FBXW7xWNT boxplot",
    "Barplots_CXCR4_pos_neg.png": "05_cxcr4_stratification.py — ranked CXCR4+ fraction barplot",
    "UMAP_CXCR4_high_low_patients.pdf": "06_dimreduction_umap.py — combined UMAP colored by patient",
    "UMAP_3Doners.pdf": "06_dimreduction_umap.py — 3-donor per-patient signature UMAP grid",
    "Palantir_CellRank_C168_trajectory.pdf": "07_trajectory_palantir.py — single-donor (C168) trajectory",
}

missing = [name for name in EXPECTED if not (fig_dir / name).exists()]
for name, desc in EXPECTED.items():
    status = "OK" if (fig_dir / name).exists() else "MISSING"
    print(f"[{status}] {name}  <-  {desc}")

if missing:
    raise SystemExit(f"\n{len(missing)} figure(s) missing: {missing}")
print(f"\nAll {len(EXPECTED)} figures present in {fig_dir}")
