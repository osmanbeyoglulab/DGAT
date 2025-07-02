import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc
import numpy as np
from scipy.spatial import Delaunay
from sklearn.metrics import f1_score,accuracy_score

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
    #print(new_adata)

    true_labels = new_adata.obs[true_label_key].values
    pred_labels = new_adata.obs[pred_label_key].values

    ari = adjusted_rand_score(true_labels, pred_labels)
    return ari    

def leiden_plot(
    adata,
    n_neighbors=10,
    resolution=0.5,
    size=5,
    points=None,
    edges=None,
    palette='tab20',
    title="Leiden Clustering by Predicted Protein",
    true_label_key='germinal_center',
    true_label_value='GC',
    None_anno_cluster=None,
    plot_type=None,
    save_fig=False
):
    """
    Leiden clustering and visualization of spatial data.
    """
    # Run Leiden clustering
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X")
    sc.tl.leiden(adata, resolution=resolution)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    if plot_type == 'F1':
        sc.pl.spatial(
            adata,
            #color=true_label_key,
            size=size,
            frameon=False,
            ax=axes[0],
            show=False,
            title="Gernimal Centers"
        )
            # Draw border lines
        if points is not None and edges is not None:
            for ii, jj in edges:
                axes[0].plot(points[[ii, jj], 0], points[[ii, jj], 1], 'k-', linewidth=1)
    else:
        if true_label_key in adata.obs:
            label_counts = adata.obs[true_label_key].value_counts()
            valid_labels = label_counts[label_counts >= 5].index
            adata_new = adata[adata.obs[true_label_key].isin(valid_labels)].copy()
        else:
            print('No label name found adata.obs')
        sc.pl.spatial(
            adata_new,
            color=true_label_key,
            size=size,
            frameon=False,
            ax=axes[0],
            show=False,
            palette = 'Paired',
            title="Pathology"
        )

    # Second subplot: Leiden + custom processing
    sc.pl.spatial(
        adata,
        color='leiden',
        size=size,
        palette=palette,
        frameon=False,
        ax=axes[1],
        show=False,
        title=title,
    )

    #  Draw border lines
    if points is not None and edges is not None:
        for ii, jj in edges:
            axes[1].plot(points[[ii, jj], 0], points[[ii, jj], 1], 'k-', linewidth=1)

    # F1 or ARI 打分
    if plot_type == 'F1':
        if true_label_key in adata.obs:
            best_score = 0
            for idx in adata.obs['leiden'].unique():
                pred = (adata.obs['leiden'] == f"{idx}").astype(int)
                true = (adata.obs[true_label_key] == true_label_value).astype(int)
                f1 = f1_score(true, pred)
                if f1 > best_score:
                    best_score = f1
            axes[1].text(
                0.02, 0.02, f"F1: {best_score:.3f}", transform=axes[1].transAxes,
                fontsize=14, color='black', ha='left', va='bottom',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7)
            )

    elif plot_type == 'ARI':
            ari_score = compute_ari_from_anndata(adata_new, true_label_key=true_label_key, pred_label_key='leiden')
            axes[1].text(
                0.02, 0.02, f"ARI: {ari_score:.3f}", transform=axes[1].transAxes,
                fontsize=20, color='black', ha='left', va='bottom',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.7)
            )

    plt.tight_layout()
    if save_fig:
        plt.savefig(f"{title.replace(' ', '_')}.png", dpi=300, bbox_inches='tight')
    plt.show()
    return adata



import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

def leiden_plot_scatter(
    adata,
    n_neighbors=10,
    resolution=0.5,
    size=100,
    true_label_key='germinal_center',
    title="Leiden Clustering by Predicted Protein",
    None_anno_cluster=None,
    save_fig=False
):
    """
    Leiden clustering and visualization of spatial data using scatter plot. This is for those samples without image related data (scalefactors and image).
    """

    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X")
    sc.tl.leiden(adata, resolution=resolution)

    coords = adata.obsm["spatial"]
    leiden_labels = adata.obs["leiden"].astype(int)
    true_labels = adata.obs[true_label_key] if true_label_key in adata.obs else None

    unique_leiden = sorted(leiden_labels.unique())
    unique_true = sorted(true_labels.unique()) if true_labels is not None else []

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    # True label
    if true_labels is not None:
        label_to_color = {
            label: plt.get_cmap("Set2")(i / max(1, len(unique_true)))
            for i, label in enumerate(unique_true)
        }
        for label in unique_true:
            mask = true_labels == label
            axes[0].scatter(coords[mask, 0], coords[mask, 1],
                            s=size, color=label_to_color[label], label=str(label), alpha=0.9)
        axes[0].legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=12, title="True Labels")
        axes[0].set_title("Pathology", fontsize=16)
    else:
        axes[0].text(0.5, 0.5, "No true labels", transform=axes[0].transAxes,
                     fontsize=14, ha='center', va='center')

    # Leiden
    for cluster in unique_leiden:
        mask = leiden_labels == cluster
        axes[1].scatter(coords[mask, 0], coords[mask, 1],
                        color=plt.get_cmap("tab20")(cluster / 20), s=size, alpha=0.9, label=f"Leiden {cluster}")
    axes[1].legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=12, title="Leiden Clusters")
    axes[1].set_title(title, fontsize=16)

    #ARI
    if true_labels is not None:
        ari_score = compute_ari_from_anndata(
            adata,
            true_label_key=true_label_key,
            pred_label_key='leiden',
            None_anno_c=None_anno_cluster
        )
        axes[1].text(
            0.02, 0.02, f"ARI: {ari_score:.3f}", transform=axes[1].transAxes,
            fontsize=20, color='black', ha='left', va='bottom',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7)
        )

    for ax in axes:
        ax.axis("equal")
        ax.axis("off")

    plt.tight_layout()
    if save_fig:
        plt.savefig(f"{title.replace(' ', '_')}_scatter.png", dpi=300, bbox_inches="tight")
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
    
    points = np.asarray(adata[adata.obs['germinal_center']=='GC'].obsm['spatial'] * 
                         adata.uns['spatial']["V1_Human_Lymph_Node"]['scalefactors']['tissue_hires_scalef'])
    points = np.vstack((points+[-r,r], points+[-r,-r], points+[r,r], points+[r,-r]))
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


