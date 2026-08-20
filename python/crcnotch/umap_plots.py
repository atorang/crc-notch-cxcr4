"""UMAP scatter-panel plotting shared by 06_dimreduction_umap.py.

Gene expression is drawn on a white->color linear colormap capped at the
90th percentile of positive cells; signature scores on a diverging
colormap capped at the 1st/99th percentile magnitude. Point size and title
format differ between the cohort-level plot (tens of thousands of cells,
s=1, title includes the vmax value) and the per-patient grid (a few
hundred to ~2000 cells per panel, s=25, title is just "{panel} ({donor})")
— both values come directly from their respective source notebooks.
"""
import numpy as np
import scipy.sparse as sp
import matplotlib.colors as mcolors


def _gene_values(adata, gene, use_raw=True):
    if use_raw:
        idx = adata.raw.var_names.get_loc(gene)
        col = adata.raw.X[:, idx]
    else:
        idx = adata.var_names.get_loc(gene)
        col = adata.X[:, idx]
    return col.toarray().flatten() if sp.issparse(col) else np.asarray(col).flatten()


def plot_gene_or_score_panel(ax, fig, adata, umap_coords, panel, is_gene, base_color,
                              point_size=1, donor=None, title_fontsize=16):
    if is_gene:
        gene_vals = _gene_values(adata, panel)
        cmap = mcolors.LinearSegmentedColormap.from_list("gene_cmap", ["#f2f1fa", base_color], N=100)
        pos_cells = gene_vals[gene_vals > 0]
        vmax = np.percentile(pos_cells, 90) if len(pos_cells) else 1
        sc_ = ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=gene_vals, cmap=cmap,
                          vmin=0, vmax=vmax, s=point_size, edgecolors="none")
        fig.colorbar(sc_, ax=ax, shrink=0.7, label="norm. expr")
        title = f"{panel} ({donor})" if donor is not None else f"{panel} (vmax 90th pct: {vmax:.2f})"
    else:
        cmap = mcolors.LinearSegmentedColormap.from_list("sig_cmap", ["#F4DF6B", "#f2f1fa", base_color], N=100)
        scores = adata.obs[panel].values.astype(float)
        vabs = max(abs(np.percentile(scores, 1)), abs(np.percentile(scores, 99)))
        sc_ = ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=scores, cmap=cmap,
                          vmin=-vabs, vmax=vabs, s=point_size, edgecolors="none")
        fig.colorbar(sc_, ax=ax, shrink=0.7, label="score")
        title = f"{panel} ({donor})" if donor is not None else panel
    ax.set_title(title, fontsize=title_fontsize)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=12 if donor is not None else 15)
    ax.set_ylabel("UMAP2", fontsize=12 if donor is not None else 15)


def plot_patient_panel(ax, fig, umap_coords, patient_labels):
    import pandas as pd
    import matplotlib.pyplot as plt

    cats = pd.Categorical(patient_labels)
    codes, uniq = cats.codes, cats.categories
    n = len(uniq)
    if n <= 10:
        cmap = plt.get_cmap("tab10", n)
    elif n <= 20:
        cmap = plt.get_cmap("tab20", n)
    else:
        colors = plt.cm.tab20b(np.linspace(0, 1, 20)).tolist() + plt.cm.tab20c(np.linspace(0, 1, 20)).tolist()
        cmap = mcolors.ListedColormap(colors[:n])
    sc_ = ax.scatter(umap_coords[:, 0], umap_coords[:, 1], c=codes, cmap=cmap,
                      vmin=-0.5, vmax=n - 0.5, s=1, edgecolors="none")
    cbar = fig.colorbar(sc_, ax=ax, shrink=0.7, ticks=range(n))
    cbar.ax.set_yticklabels(list(uniq), fontsize=6)
    cbar.set_label("Patient ID")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=12); ax.set_ylabel("UMAP2", fontsize=12)
    ax.set_title("donor_id", fontsize=16)
