"""07_trajectory_palantir.py

Palantir trajectory inference + CellRank pseudotime projection for donor
C168.

Root cell: closest-to-centroid of the top 20 cells by Fetal_curated score.
Terminal cell: closest-to-centroid of the top 20 cells by LGR5 expression.

Input:  data/processed/pelka/Pelka_CXCR4+_processed_C168.h5ad (written by 06)
Output: results/figures/Palantir_CellRank_C168_trajectory.pdf

"""
import io
import sys
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "python"))
from crcnotch.trajectory import run_cellrank_pseudotime_kernel, run_palantir_trajectory  # noqa: E402

GENE_PANEL = ["CXCR4", "LGR5", "AXIN2", "NKD1", "ANXA3", "ANXA1", "TACSTD2"]
PANEL_CMAP = "Spectral_r"


def plot_gene_expression_premium(df_expr, branch_name, color_trop2="#c02ded", color_lgr5="#2D5A43",
                                  color_cxcr4="#08A047", linewidth=1, figsize=(4.4, 3.4)):
    """Min-max-normalized MAGIC gene trends for TACSTD2/LGR5/CXCR4 along one
    branch — matches the original notebook's plot_gene_expression_premium()."""
    df_branch = df_expr[branch_name] if branch_name in df_expr.columns.get_level_values(0) else df_expr
    pseudotime = df_branch.index.values.astype(float)

    def get_gene(df, name):
        match = [c for c in df.columns if str(c).lower() == name.lower()]
        if not match:
            raise KeyError(f"Gene '{name}' not found. Available: {list(df.columns)}")
        return df[match[0]].values.flatten()

    def min_max_scale(x):
        xmin, xmax = np.nanmin(x), np.nanmax(x)
        return (x - xmin) / (xmax - xmin) if xmax > xmin else np.zeros_like(x)

    tacstd2 = min_max_scale(get_gene(df_branch, "Tacstd2"))
    lgr5 = min_max_scale(get_gene(df_branch, "Lgr5"))
    cxcr4 = min_max_scale(get_gene(df_branch, "Cxcr4"))

    fig, ax = plt.subplots(figsize=figsize, dpi=200)
    for expr, label, color in [(tacstd2, "Tacstd2", color_trop2), (lgr5, "Lgr5", color_lgr5),
                                (cxcr4, "Cxcr4", color_cxcr4)]:
        ax.plot(pseudotime, expr, label=label, color=color, linewidth=linewidth, zorder=3)
        ax.fill_between(pseudotime, expr, color=color, alpha=0.08, zorder=2)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#cccccc")
    ax.spines["bottom"].set_color("#cccccc")
    ax.grid(axis="y", linestyle="--", alpha=0.4, color="#cccccc", zorder=1)
    ax.set_xlabel("Trajectory Pseudotime", fontsize=10, labelpad=6, color="#333333")
    ax.set_ylabel("Relative Expression", fontsize=10, labelpad=6, color="#333333")
    ax.tick_params(colors="#555555", labelsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.02, 1.05)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=3, frameon=False,
              fontsize=9, handletextpad=0.5, columnspacing=1.5)
    fig.tight_layout()
    return fig


def fig_to_array(fig, dpi=200):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    buf.seek(0)
    return plt.imread(buf)


def render_trajectory_panel(adata, donor, figsize=(4.4, 4.0)):
    import palantir as pt

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    pt.plot.plot_trajectories(adata, pseudotime_interval=(0, 1.0), lw=0.5,
                               scanpy_kwargs={"size": 25, "cmap": PANEL_CMAP}, ax=ax)
    if len(fig.axes) > 1:
        cax = fig.axes[-1]
        pos = cax.get_position()
        new_height = pos.height * 0.5
        new_y0 = pos.y0 + (pos.height - new_height) / 2
        cax.set_position([pos.x0, new_y0, pos.width, new_height])
    ax.set_title(f"Patient {donor}", fontsize=12)
    return fig


def render_heatmap_panel(adata, width=4.4, row_height=0.4):
    import palantir as pt

    return pt.plot.plot_gene_trend_heatmaps(adata, GENE_PANEL, cmap="coolwarm", basefigsize=(width, row_height))


def render_expression_panel(adata, figsize=(4.4, 3.4)):
    branches = list(adata.obsm["palantir_fate_probabilities"].columns)
    n_points = adata.varm[f"gene_trends_{branches[0]}"].shape[1]
    pseudotime_grid = np.linspace(0, 1, n_points)
    df_expr = pd.concat({br: adata.varm[f"gene_trends_{br}"].loc[GENE_PANEL].T for br in branches}, axis=1)
    df_expr.index = pseudotime_grid
    df_expr.index.name = "pseudotime"
    return plot_gene_expression_premium(df_expr, branch_name=branches[0], figsize=figsize)


def render_cellrank_panel(kernel, figsize=(4.4, 4.0)):
    fig, ax = plt.subplots(figsize=figsize)
    kernel.plot_projection(basis="X_umap", color="palantir_pseudotime", color_map=PANEL_CMAP,
                            size=25, legend_loc=None, colorbar=True, alpha=1.0, ax=ax, show=False)
    return fig


def main():
    cfg = yaml.safe_load((REPO_ROOT / "config" / "config.yaml").read_text())
    traj = cfg["trajectory_palantir"]
    seed = cfg["seed"]

    pelka_dir = REPO_ROOT / cfg["paths"]["data_processed"] / "pelka"
    fig_dir = REPO_ROOT / cfg["paths"]["results_figures"]
    fig_dir.mkdir(parents=True, exist_ok=True)

    donor = traj["donor"]
    adata = sc.read_h5ad(pelka_dir / f"Pelka_CXCR4+_processed_{donor}.h5ad")
    adata.X = adata.raw.X.copy()
    expr = adata[:, "LGR5"].X
    adata.obs["LGR5"] = (expr.toarray() if hasattr(expr, "toarray") else np.asarray(expr)).flatten()

    adata, root_cell, terminal_cells = run_palantir_trajectory(
        adata, root_signature=traj["root_signature"], terminal_signature=traj["terminal_signature"],
        n_components=traj["n_diffusion_components"], knn=traj["knn"], num_waypoints=traj["num_waypoints"],
        root_top_n=traj["root_top_n_cells"], terminal_top_n=20, terminal_n_clusters=1, seed=seed,
    )
    print(f"{donor}: n_cells={adata.n_obs}, root={root_cell}, terminal={terminal_cells}")
    print(f"pseudotime range: {adata.obs['palantir_pseudotime'].min()} - {adata.obs['palantir_pseudotime'].max()}")

    kernel = run_cellrank_pseudotime_kernel(adata, time_key="palantir_pseudotime")

    import palantir as pt

    # Required before compute_gene_trends: assigns adata.obsm["branch_masks"].
    masks = pt.presults.select_branch_cells(adata)
    if not isinstance(masks, pd.DataFrame):
        masks = pd.DataFrame(masks, index=adata.obs_names, columns=adata.obsm["palantir_fate_probabilities"].columns)
    adata.uns["palantir_branch_masks"] = masks

    # Rename the terminal-state branch 
    old_branches = list(adata.obsm["palantir_fate_probabilities"].columns)
    branch_map = {old: letter for old, letter in zip(old_branches, "ABCDEF")}
    adata.obsm["palantir_fate_probabilities"].rename(columns=branch_map, inplace=True)
    adata.uns["palantir_branch_masks"].rename(columns=branch_map, inplace=True)
    adata.obsm["branch_masks"].rename(columns=branch_map, inplace=True)

    pt.presults.compute_gene_trends(adata, expression_key="MAGIC_imputed_data", pseudo_time_key="palantir_pseudotime")

    row_builders = [
        ("Palantir trajectory", lambda: render_trajectory_panel(adata, donor)),
        ("Gene-trend heatmap", lambda: render_heatmap_panel(adata)),
        ("Branch expression", lambda: render_expression_panel(adata)),
        ("CellRank projection", lambda: render_cellrank_panel(kernel)),
    ]
    panel_images = {label: fig_to_array(builder()) for label, builder in row_builders}

    col_width = 4.4
    row_labels = [label for label, _ in row_builders]
    row_heights = [col_width * (panel_images[label].shape[0] / panel_images[label].shape[1]) for label in row_labels]

    fig = plt.figure(figsize=(col_width + 0.9, sum(row_heights) + 0.6))
    gs = gridspec.GridSpec(len(row_labels), 1, figure=fig, height_ratios=row_heights,
                            left=0.09, right=0.99, top=0.94, bottom=0.02, hspace=0.15)
    for r, label in enumerate(row_labels):
        ax = fig.add_subplot(gs[r, 0])
        ax.imshow(panel_images[label])
        ax.set_xticks([]); ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_ylabel(label, fontsize=11, rotation=90, labelpad=8)

    fig.suptitle(f"Palantir / CellRank trajectory {donor}", fontsize=14)
    out_path = fig_dir / "Palantir_CellRank_C168_trajectory.pdf"
    fig.savefig(out_path, format="pdf")
    plt.close(fig)
    print("wrote", out_path)


if __name__ == "__main__":
    main()
