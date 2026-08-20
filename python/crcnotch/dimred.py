"""HVG -> scale -> PCA -> neighbors -> Leiden -> UMAP pipeline, shared by
the cohort-level (config: dimensionality_reduction.cohort_level) and
per-patient (config: dimensionality_reduction.per_patient) stages. The
log-normalized layer is preserved in `.raw` before scaling, since
signature_scoring.score_signature_zscore needs it afterward.
"""


def run_dimred(adata, n_hvg, n_pcs, n_neighbors, leiden_resolution, seed,
               hvg_flavor="seurat", scale_clip_max=10, leiden_flavor="igraph",
               leiden_n_iterations=2, umap_min_dist=None, pca_solver="arpack",
               n_pcs_for_neighbors=None):
    """`n_pcs` PCA components are computed; by default all of them feed the
    neighbor graph. Pass `n_pcs_for_neighbors` < `n_pcs` to compute more
    components than are actually used for the graph — the per-patient step
    does this (20 computed, 10 used), matching the original notebook."""
    import scanpy as sc

    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor=hvg_flavor)
    adata.raw = adata.copy()  # preserve log-normalized layer for signature scoring
    sc.pp.scale(adata, max_value=scale_clip_max)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver=pca_solver, random_state=seed)
    neighbor_pcs = n_pcs_for_neighbors if n_pcs_for_neighbors is not None else n_pcs
    sc.pp.neighbors(adata, n_pcs=neighbor_pcs, n_neighbors=n_neighbors, random_state=seed)
    sc.tl.leiden(adata, resolution=leiden_resolution, flavor=leiden_flavor,
                 n_iterations=leiden_n_iterations, random_state=seed)
    umap_kwargs = {"random_state": seed}
    if umap_min_dist is not None:
        umap_kwargs["min_dist"] = umap_min_dist
    sc.tl.umap(adata, **umap_kwargs)
    return adata
