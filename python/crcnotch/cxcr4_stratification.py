"""Per-patient CXCR4 stratification: fraction of CXCR4-positive malignant
cells per patient, high/low/intermediate grouping, and the ranked-fraction
barplot. 
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp


def compute_cxcr4_positivity(adata, gene_symbol_col="gene_symbol", target_gene="CXCR4",
                              positivity_threshold=0.0):
    """Flag each cell CXCR4+ if its log-normalized expression exceeds
    `positivity_threshold` (config: cxcr4_stratification.positivity_threshold).
    """
    gene_ids = adata.var_names[adata.var[gene_symbol_col] == target_gene].tolist()
    if not gene_ids:
        raise ValueError(f"{target_gene} not found in var[{gene_symbol_col}]")
    gene_id = gene_ids[0]
    expr = adata[:, gene_id].X
    expr = np.asarray(expr.todense()).ravel() if sp.issparse(expr) else np.asarray(expr).ravel()
    adata.obs["CXCR4_expr"] = expr
    adata.obs[f"{target_gene}_pos"] = expr > positivity_threshold
    return adata


def per_patient_fraction(adata, patient_key, min_tumor_cells, positivity_col="CXCR4_pos",
                          expr_col="CXCR4_expr"):
    per_pat = (
        adata.obs.groupby(patient_key)
        .agg(n_tumor=(positivity_col, "size"),
             n_cxcr4_pos=(positivity_col, "sum"),
             mean_cxcr4=(expr_col, "mean"))
        .assign(frac_cxcr4_pos=lambda d: d["n_cxcr4_pos"] / d["n_tumor"])
    )
    per_pat = per_pat[per_pat["n_tumor"] >= min_tumor_cells].sort_values("frac_cxcr4_pos", ascending=False)
    return per_pat


def stratify_patients(per_pat, high_cutoff, low_cutoff, frac_col="frac_cxcr4_pos"):
    conditions = [per_pat[frac_col] >= high_cutoff, per_pat[frac_col] <= low_cutoff]
    choices = ["CXCR4_high", "CXCR4_low"]
    per_pat = per_pat.copy()
    per_pat["CXCR4_group"] = np.select(conditions, choices, default="CXCR4_med")
    return per_pat


def plot_ranked_fraction_barplot(per_pat, high_cutoff, low_cutoff, out_path, frac_col="frac_cxcr4_pos"):
    import matplotlib.pyplot as plt

    group_colors = {"CXCR4_high": "#c0392b", "CXCR4_med": "#bdc3c7", "CXCR4_low": "#2980b9"}
    colors = per_pat["CXCR4_group"].map(group_colors)

    fig, ax = plt.subplots(figsize=(11, 4))
    ax.bar(range(len(per_pat)), per_pat[frac_col], color=colors)
    ax.axhline(high_cutoff, ls="--", c="#c0392b", lw=1.2, label=f"High Threshold (>= {high_cutoff:.2f})")
    ax.axhline(low_cutoff, ls="--", c="#2980b9", lw=1.2, label=f"low Threshold (<= {low_cutoff:.2f})")
    ax.set_xticks(range(len(per_pat)))
    ax.set_xticklabels(per_pat.index, rotation=90, fontsize=6)
    ax.set_ylabel("fraction CXCR4+ tumor cells")
    ax.set_title("Per-patient CXCR4+ tumor-cell fraction")
    ax.legend()
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path
