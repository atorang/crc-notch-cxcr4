"""QC filtering shared by every scRNA-seq preprocessing script.

Reproduces the original notebooks' QC block exactly: cells with fewer than
`min_genes` or more than `max_genes` detected genes are dropped (the upper
bound is a doublet/multiplet heuristic), cells with more than `max_pct_mt`
mitochondrial-gene counts are dropped, then genes detected in fewer than
`min_cells_per_gene` cells are dropped. Metagene "pseudo-genes" introduced
by the Pelka et al. paper (`pEpi*` columns, matched case-insensitively as
`^p[A-Z]`) are extracted into `.obs` before being removed from `.var`, since
several downstream notebooks read them as per-cell metaprogram scores.
"""
import re

import numpy as np
import pandas as pd
import scipy.sparse as sp


def extract_pelka_metaprograms(adata):
    """Move Pelka et al.'s `pEpi*`-style metagene pseudo-genes from `.var`
    into `.obs`, then drop them from the expression matrix so they don't
    contaminate HVG selection or signature scoring.
    """
    is_metagene = adata.var_names.str.match(r"^p[A-Z]")
    target_states = adata.var_names[is_metagene]
    new_cols = {}
    for state in target_states:
        idx = adata.var_names.get_loc(state)
        col = adata.X[:, idx]
        new_cols[state] = col.toarray().flatten() if sp.issparse(col) else np.asarray(col).flatten()
    adata.obs = pd.concat([adata.obs, pd.DataFrame(new_cols, index=adata.obs_names)], axis=1)
    return adata[:, ~is_metagene].copy()


def qc_filter(adata, min_genes, max_genes, max_pct_mt, min_cells_per_gene, mt_prefix="MT-"):
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    import scanpy as sc
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    valid_genes = (adata.obs["n_genes_by_counts"] >= min_genes) & (adata.obs["n_genes_by_counts"] <= max_genes)
    valid_mt = adata.obs["pct_counts_mt"] <= max_pct_mt
    adata = adata[valid_genes & valid_mt, :].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    return adata


def ensure_log_normalized(adata, raw_max_threshold=30):
    """Log1p-transform `adata.X` only if it still looks count-like (max
    value above `raw_max_threshold`); otherwise assume it is already
    log-normalized. """
    import scanpy as sc
    X = adata.X
    xmax = X.max() if not sp.issparse(X) else X.max()
    if xmax > raw_max_threshold:
        sc.pp.log1p(adata)
    return adata
