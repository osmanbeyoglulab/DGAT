import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc
import numpy as np
from scipy.spatial import Delaunay
from sklearn.metrics import f1_score, accuracy_score

from sklearn.metrics import adjusted_rand_score


def compute_ari_from_anndata(adata, true_label_key='Pathology', pred_label_key='leiden', None_anno_c=None):
    label_counts = adata.obs[true_label_key].value_counts()
    valid_labels = label_counts[label_counts >= 5].index
    adata = adata[adata.obs[true_label_key].isin(valid_labels)].copy()

    if None_anno_c:
        new_adata = adata[adata.obs[true_label_key] != None_anno_c, :].copy()
    else:
        new_adata = adata.copy()

    mask = (
            new_adata.obs[true_label_key].notna() &
            new_adata.obs[pred_label_key].notna()
    )
    new_adata = new_adata[mask, :].copy()
    # print(new_adata)

    true_labels = new_adata.obs[true_label_key].values
    pred_labels = new_adata.obs[pred_label_key].values

    ari = adjusted_rand_score(true_labels, pred_labels)
    return ari


def leiden_plot(adata, n_neighbors=10, resolution=0.5, size=5, points=None, edges=None, palette='tab20',
                title="Leiden Clustering by Predicted Protein", true_label_key='germinal_center', true_label_value='GC',
                None_anno_cluster=None, plot_type=None, save_fig=False):
    # Run Leiden clustering
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X")
    sc.tl.leiden(adata, resolution=resolution)
    fig, ax = plt.subplots(figsize=(5, 5))

    sc.pl.spatial(
        adata,
        color='leiden',
        size=size,
        palette=palette,
        frameon=False,
        ax=ax,
        show=False,
        # title=method_list[idx],
        # legend_loc=None,
    )
    if points is not None and edges is not None:
        for ii, jj in edges:
            ax.plot(points[[ii, jj], 0], points[[ii, jj], 1], 'k-', linewidth=1)
    if plot_type == 'F1':
        if true_label_key in adata.obs:
            best_score = 0
            for idx in adata.obs['leiden'].unique():
                pred = (adata.obs['leiden'] == f"{idx}").astype(int)
                true = (adata.obs[true_label_key] == true_label_value).astype(int)
                f1 = f1_score(true, pred)
                if f1 > best_score:
                    best_score = f1
            ax.text(
                0.02, 0.02, f"F1: {best_score:.3f}", transform=ax.transAxes,
                fontsize=14, color='black', ha='left', va='bottom',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7)
            )
    elif plot_type == 'ARI':
        if true_label_key in adata.obs:
            label_counts = adata.obs[true_label_key].value_counts()
            valid_labels = label_counts[label_counts >= 5].index
            # remove low number case
            adata_new = adata[adata.obs[true_label_key].isin(valid_labels)].copy()
            print(adata_new)
            ari_score = compute_ari_from_anndata(adata_new, true_label_key=true_label_key, pred_label_key='leiden')
            ax.text(
                0.02, 0.02, f"ARI: {ari_score:.3f}", transform=ax.transAxes,
                fontsize=20, color='black', ha='left', va='bottom',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7)
            )
            # print(ari_score)
    # plt.close()
    plt.show()
    return adata


def leiden_plot_scatter(adata, n_neighbors=10, resolution=0.5, size=100, true_label_key='germinal_center',
                        title="Leiden Clustering by Predicted Protein", None_anno_cluster=None, save_fig=False):
    # Run Leiden clustering
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X")
    sc.tl.leiden(adata, resolution=resolution)

    # Ensure spatial coordinates are correctly set
    if "X_spatial" not in adata.obsm:
        adata.obsm["X_spatial"] = adata.obsm["spatial"]

    # Extract spatial coordinates and Leiden clusters
    coords = adata.obsm["X_spatial"]
    labels = adata.obs["leiden"].astype(int)
    unique_labels = sorted(labels.unique())

    # Create the scatter plot
    plt.figure(figsize=(8, 8))
    for cluster in unique_labels:
        mask = labels == cluster
        plt.scatter(coords[mask, 0], coords[mask, 1], color=plt.get_cmap("tab20")(cluster / 20), s=size, alpha=0.9,
                    label=f"Leiden {cluster}")

    # Add title and axis formatting
    plt.title(title, fontsize=16)
    plt.axis("equal")
    plt.xticks([])
    plt.yticks([])
    plt.axis("off")
    # Create a custom legend
    handles = [mpatches.Patch(color=plt.get_cmap("tab20")(cluster / 20), label=f"Leiden {cluster}") for cluster in
               unique_labels]
    plt.legend(handles=handles, loc="center left", bbox_to_anchor=(1, 0.5), fontsize=12, title="Leiden Clusters")
    if true_label_key in adata.obs:
        ari_score = compute_ari_from_anndata(adata, true_label_key=true_label_key, pred_label_key='leiden',
                                             None_anno_c=None_anno_cluster)
        plt.text(
            0.02, 0.02, f"ARI: {ari_score:.3f}",
            fontsize=20, color='black', ha='left', va='bottom',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7)
        )
    # Save if requested
    if save_fig:
        plt.savefig(f"{title.replace(' ', '_')}.png", dpi=300, bbox_inches="tight")

    plt.show()
    return adata


def adata_to_df(adata):
    if isinstance(adata.X, np.ndarray):
        X = adata.X
    else:
        X = adata.X.toarray()
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    return df


def find_edges_LN(adata, r=10, alpha=20, only_outer=True):
    points = np.asarray(adata[adata.obs['germinal_center'] == 'GC'].obsm['spatial'] *
                        adata.uns['spatial']["V1_Human_Lymph_Node"]['scalefactors']['tissue_hires_scalef'])
    points = np.vstack((points + [-r, r], points + [-r, -r], points + [r, r], points + [r, -r]))
    assert points.shape[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        if (i, j) in edges or (j, i) in edges:
            assert (j, i) in edges, "Can't go twice over same directed edge right?"
            if only_outer:
                edges.remove((j, i))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    edges = set()
    for ia, ib, ic in tri.simplices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        a = np.sqrt((pa[0] - pb[0]) ** 2 + (pa[1] - pb[1]) ** 2)
        b = np.sqrt((pb[0] - pc[0]) ** 2 + (pb[1] - pc[1]) ** 2)
        c = np.sqrt((pc[0] - pa[0]) ** 2 + (pc[1] - pa[1]) ** 2)
        s = (a + b + c) / 2.0
        area = np.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)
    return points, edges


