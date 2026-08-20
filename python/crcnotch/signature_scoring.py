from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def score_signature_zscore(adata, gene_sets, use_raw=True, layer=None, verbose=True):
    """Score each signature in `gene_sets` on every cell in `adata`.

    For each signature: take the genes present in the scoring layer, compute
    each gene's Z-score across the cells being scored, then average the per-gene Z-scores for
    each cell. Writes one column per signature into `adata.obs`.

    Parameters
    ----------
    adata : AnnData
        The population of cells to score. Z-scoring is computed across
        exactly the cells in this object, pass a per-patient subset for a
        within-patient score, or the full cohort for a global score.
    gene_sets : dict[str, list[str]]
        Signature name -> list of gene symbols.
    use_raw : bool
        If True (default) and `adata.raw` is set, score on `adata.raw.X`
        (the log-normalized reference layer). If `adata.raw` is not set, falls back to `adata.X` with a
        warning.
    layer : str or None
        If given, score on `adata.layers[layer]` instead of raw/`.X`.
    """
    if layer is not None:
        var_names = adata.var_names
        X = adata.layers[layer]
    elif use_raw and adata.raw is not None:
        var_names = adata.raw.var_names
        X = adata.raw.X
    else:
        if use_raw:
            print("[score_signature_zscore] WARNING: use_raw=True but adata.raw is not set, "
                  "scoring on adata.X, which may be the scaled/clipped layer, not log-normalized.")
        var_names = adata.var_names
        X = adata.X

    for name, genes in gene_sets.items():
        genes_found = [g for g in genes if g in var_names]
        genes_missing = len(genes) - len(genes_found)
        if not genes_found:
            if verbose:
                print(f"[{name}] No genes found, skipping")
            continue
        if verbose:
            print(f"[{name}] {len(genes_found)}/{len(genes)} genes found ({genes_missing} missing)")
        idx = [var_names.get_loc(g) for g in genes_found]
        X_sub = X[:, idx]
        if sp.issparse(X_sub):
            X_sub = X_sub.toarray()
        X_sub = np.asarray(X_sub)
        means = X_sub.mean(axis=0)
        stds = X_sub.std(axis=0)
        stds[stds == 0] = 1
        X_z = (X_sub - means) / stds
        adata.obs[name] = X_z.mean(axis=1)
    if verbose:
        print("Done! Scores saved to adata.obs")


def load_gene_sets(cfg, repo_root):
    """Build the gene_sets dict from config.yaml's signature_scoring.gene_sets
    entry: inline lists (Notch_Sander, Fetal_curated, WNT_curated) are used; 
    entries with a `file` key are read relative to `repo_root`
    (e.g. "data/metadata/signatures/fetal_mustata.tsv").
    """
    repo_root = Path(repo_root)
    gene_sets = {}
    for name, spec in cfg["signature_scoring"]["gene_sets"].items():
        if "genes" in spec:
            gene_sets[name] = list(spec["genes"])
        elif "file" in spec:
            df = pd.read_csv(repo_root / spec["file"], sep="\t")
            gene_sets[name] = df["gene"].astype(str).tolist()
        else:
            raise ValueError(f"gene set '{name}' has neither 'genes' nor 'file' in config.yaml")
    return gene_sets
