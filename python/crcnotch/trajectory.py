"""Palantir trajectory inference + CellRank pseudotime projection for a
single donor
"""
import numpy as np
import pandas as pd


def get_root_cell(adata, cell_subset, coord_key="DM_EigenVectors_multiscaled"):
    """Closest cell to the centroid of `cell_subset` in Palantir's
    multiscale diffusion space."""
    dm_coords = adata[cell_subset].obsm[coord_key]
    centroid = dm_coords.mean(axis=0)
    dists = np.linalg.norm(dm_coords - centroid, axis=1)
    return cell_subset[np.argmin(dists)]


def get_terminal_cells(adata, signature_key, top_n, n_clusters, seed,
                        coord_key="DM_EigenVectors_multiscaled"):
    """Cluster the top-`top_n` cells by `signature_key` into `n_clusters`
    groups in diffusion space and return each cluster's medoid as a
    terminal-state candidate."""
    from sklearn.cluster import KMeans

    top_cells = adata.obs.nlargest(top_n, signature_key).index
    dm_coords = adata[top_cells].obsm[coord_key]
    labels = KMeans(n_clusters=n_clusters, random_state=seed).fit(dm_coords).labels_

    terminal_cells = []
    for cluster_id in range(n_clusters):
        mask = labels == cluster_id
        cluster_cells = top_cells[mask]
        cluster_coords = dm_coords[mask]
        centroid = cluster_coords.mean(axis=0)
        dists = np.linalg.norm(cluster_coords - centroid, axis=1)
        terminal_cells.append(cluster_cells[np.argmin(dists)])
    return terminal_cells


def run_palantir_trajectory(adata, root_signature, terminal_signature, n_components, knn,
                             num_waypoints, root_top_n, terminal_top_n, terminal_n_clusters, seed):
    import palantir as pt

    pt.utils.run_diffusion_maps(adata, pca_key="X_pca", n_components=n_components, knn=knn, seed=seed)
    pt.utils.determine_multiscale_space(adata)
    pt.utils.run_magic_imputation(adata)

    root_candidates = adata.obs.nlargest(root_top_n, root_signature).index
    root_cell = get_root_cell(adata, root_candidates)
    terminal_cells = get_terminal_cells(adata, terminal_signature, terminal_top_n, terminal_n_clusters, seed)

    pt.core.run_palantir(
        adata, early_cell=root_cell, terminal_states=terminal_cells,
        num_waypoints=num_waypoints, use_early_cell_as_start=True, knn=knn, seed=seed,
    )
    return adata, root_cell, terminal_cells


def run_cellrank_pseudotime_kernel(adata, time_key="palantir_pseudotime"):
    import cellrank as cr

    pk = cr.kernels.PseudotimeKernel(adata, time_key=time_key)
    pk.compute_transition_matrix(threshold_scheme="soft")
    return pk
