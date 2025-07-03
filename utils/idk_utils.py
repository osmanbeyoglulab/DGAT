import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import scanpy as sc
import numpy as np
from scipy.spatial import Delaunay
from sklearn.metrics import f1_score,accuracy_score
from sklearn.metrics import adjusted_rand_score
from statsmodels.stats.multitest import multipletests
from scipy.stats import t
import seaborn as sns
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

# leiden plot for predicted protein with true protein

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

def find_edges(adata, r=10, alpha=15, only_outer=True, cluster_name = 'germinal_center', key_name = 'GC'):
    points = np.asarray(adata[adata.obs[cluster_name] == key_name].obsm['spatial'] *
                        adata.uns['spatial'][list(adata.uns['spatial'].keys())[0]]['scalefactors'][
                            'tissue_hires_scalef'])
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

# Evaluation by mRNA, Measured Protein, and Predicted Protein
def leiden_plot_eva(
        adata,
        pdata,
        predata,
        points=None,
        edges=None,
        n_neighbor_list=None,
        resolution_list=None,
        size=5,
        palette='tab20',
        title="Leiden Clustering by Predicted Protein",
        true_label_key='germinal_center',
        true_label_value='GC',
        None_anno_cluster=None,
        title_list = ['mRNA','Measured Protein', 'Predicted Protein'],
        plot_type=None,

        save_fig=False

):
    """
    Leiden clustering and visualization of spatial data.
    """
    if n_neighbor_list is None:
        n_neighbor_list = [10, 10, 10]
    if resolution_list is None:
        resolution_list = [0.3, 0.5, 0.7]

    data_list = [adata,pdata,predata]
    for data, n_neighbors, resolution in zip(data_list,n_neighbor_list,resolution_list):
        sc.pp.neighbors(data, n_neighbors=n_neighbors, use_rep="X")
        sc.tl.leiden(data, resolution=resolution)



    # First column: True label
    if plot_type == 'F1':
        fig, axes = plt.subplots(1, 4, figsize=(20, 6))
        sc.pl.spatial(
            adata,
            size=size,
            frameon=False,
            ax=axes[0],
            show=False,
            title="Gernimal Centers"
        )
        if points is not None and edges is not None:
            for ii, jj in edges:
                axes[0].plot(points[[ii, jj], 0], points[[ii, jj], 1], 'k-', linewidth=1)

    #最后三个axes:
    for i,ax in enumerate(axes[-3:]):
        ax.axis("off")
        sc.pl.spatial(
            data_list[i],
            color='leiden',
            size=size,
            frameon=False,
            ax=ax,
            show=False,
            palette='Paired',
            title=title_list[i]
        )
        if points is not None and edges is not None:
            for ii, jj in edges:
                ax.plot(points[[ii, jj], 0], points[[ii, jj], 1], 'k-', linewidth=1)

        if plot_type=='F1':
            if true_label_key in data_list[i].obs:
                best_score = 0
                for idx in data_list[i].obs['leiden'].unique():
                    pred = (data_list[i].obs['leiden'] == f"{idx}").astype(int)
                    true = (data_list[i].obs[true_label_key] == true_label_value).astype(int)
                    f1 = f1_score(true, pred)
                    if f1 > best_score:
                        best_score = f1
                ax.text(
                    0.02, 0.02, f"F1: {best_score:.3f}", transform=ax.transAxes,
                    fontsize=14, color='black', ha='left', va='bottom',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7)
                )


    if save_fig:
        plt.savefig(f"{title.replace(' ', '_')}.png", dpi=300, bbox_inches='tight')
    plt.tight_layout()
    plt.show()



def plot_spatial_expression(ad, genes, points, edges, method, figsize= 2.5,fontsize = 9):
    df = sc.get.rank_genes_groups_df(ad, group='GC')
    df.index = df['names']

    ngenes = len(genes)
    fig, axs = plt.subplots(1, ngenes, figsize=(ngenes*figsize, figsize), dpi=150)
    plt.rc('font', size=fontsize)
    for i in range(ngenes):
        gene = genes[i]
        sc.pl.spatial(ad, color=gene,
            size=1.8, alpha_img=0, color_map="plasma", ax=axs[i], show=False,
            legend_fontsize=fontsize, colorbar_loc='right')
        title = f'{method}-{gene}\np_adj=%.2e'%df.loc[gene,'pvals_adj']
        axs[i].set_title(title, fontsize=fontsize)
        axs[i].set_xlabel("")
        axs[i].set_ylabel("")
        for ii, jj in edges:
            axs[i].plot(points[[ii, jj], 0], points[[ii, jj], 1], 'k', linewidth=0.8)
    plt.tight_layout(pad=0.6)


def infer_celltype_activity(adata):
    A = adata.obsm['celltype_major'].to_numpy()
    b = adata.to_df().to_numpy()
    cov = np.dot(b.T - b.mean(), A - A.mean(axis=0)) / (b.shape[0] - 1)
    ssd = np.std(A, axis=0, ddof=1) * np.std(b, axis=0, ddof=1).reshape(-1, 1)
    r = cov / ssd

    n_samples = b.shape[1]
    n_features, n_fsets = A.shape
    df = n_features - 2
    es = r * np.sqrt(df / ((1.0 - r + 1.0e-16) * (1.0 + r + 1.0e-16)))

    pv = t.sf(abs(es), df) * 2

    cts = adata.obsm['celltype_major'].columns
    tfs = adata.to_df().columns

    estimate = pd.DataFrame(es, index=tfs, columns=cts)
    estimate.name = 'coef'
    pvals = pd.DataFrame(pv, index=tfs, columns=cts)
    pvals.name = 'pvals'

    df_1 = estimate.stack().reset_index()
    df_1.columns = ['pt', 'ct', 'coef']

    df_2 = pvals.stack().reset_index()
    df_2.columns = ['pt', 'ct', 'pval']

    df_ct_tf = pd.merge(df_1, df_2, on=['pt', 'ct'])
    df_ct_tf["p_adj"] = multipletests(df_ct_tf['pval'], alpha=0.01, method="fdr_bh")[1]
    df_ct_tf["neg_log_p_adj"] = -np.log10(df_ct_tf["p_adj"] + 1e-100)
    return df_ct_tf

def merge_celltypes(adata):
    df_celltype = pd.DataFrame(0, index=adata.obsm['celltype'].index,
                          columns = ['B_Cycling', 'B_GC', 'B_IFN', 'B_activated', 'B_mem', 'B_naive', 'B_plasma', 'B_preGC',
                                     'DC', 'Endo', 'FDC', 'ILC', 'Macrophages', 'Mast', 'Monocytes', 'NK', 'NKT',
                                     'T_CD4+', 'T_CD8+', 'T_Treg', 'T_TIM3+', 'T_TfR', 'VSMC'])
    for ct in ['B_Cycling', 'B_IFN', 'B_activated', 'B_mem', 'B_naive', 'B_plasma', 'B_preGC',
               'Endo', 'FDC', 'ILC', 'Mast', 'Monocytes', 'NK', 'NKT', 'T_Treg', 'T_TIM3+', 'T_TfR', 'VSMC']:
        df_celltype[ct] = adata.obsm['celltype'][ct]

    for ct in ['B_GC_DZ', 'B_GC_LZ', 'B_GC_prePB']:
        df_celltype['B_GC'] += adata.obsm['celltype'][ct]

    for ct in ['DC_CCR7+', 'DC_cDC1', 'DC_cDC2', 'DC_pDC']:
        df_celltype['DC'] += adata.obsm['celltype'][ct]

    for ct in ['Macrophages_M1', 'Macrophages_M2']:
        df_celltype['Macrophages'] += adata.obsm['celltype'][ct]

    for ct in ['T_CD4+', 'T_CD4+_TfH', 'T_CD4+_TfH_GC', 'T_CD4+_naive']:
        df_celltype['T_CD4+'] += adata.obsm['celltype'][ct]

    for ct in ['T_CD8+_CD161+', 'T_CD8+_cytotoxic', 'T_CD8+_naive']:
        df_celltype['T_CD8+'] += adata.obsm['celltype'][ct]
    return df_celltype

def plot_heatmap(df_ct_tf,  pt_list, ct_list,clip=10):
    data = df_ct_tf.query("pt in @pt_list and ct in @ct_list")
    data.columns = ['pt', 'ct', 'Cell type-specific protein Score', 'pval', 'p_adj', '-log(p_adj)']
    x = 'pt'
    y = 'ct'
    hue = 'Cell type-specific protein Score'

    data[x] = data[x].astype("category")
    data[y] = data[y].astype("category")
    x_lab = data[x].cat.categories
    y_lab = data[y].cat.categories

    f = sns.clustermap(data.pivot(index=y, columns=x, values=hue),figsize=(0.1,0.1), cmap='PiYG')
    x_lab = x_lab[f.dendrogram_col.reordered_ind]
    y_lab = y_lab[f.dendrogram_row.reordered_ind]
    # print(x_lab)
    # print(y_lab)

    data[x] = data[x].cat.reorder_categories(x_lab)
    data[y] = data[y].cat.reorder_categories(y_lab)
    data = data.sort_values([x, y])
    data[hue] = data[hue].clip(-clip, clip)

    figsize = 0.2
    plt.figure(figsize=(figsize*len(x_lab), figsize*len(y_lab)), dpi=150)
    plt.rc('font', size=9)
    ax = sns.scatterplot(data=data,x=x, y=y, palette="PiYG_r", hue=hue, size="-log(p_adj)")
    plt.legend(bbox_to_anchor=(1.5,1.), loc='upper right',
        columnspacing=0.5, handletextpad=0, frameon=False, fontsize=9)

    ax.set_xticklabels(x_lab, rotation = 90)
    ax.set_xlim(-0.5, -0.5+len(x_lab))
    ax.set_ylim(-0.5, -0.5+len(y_lab))
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.close(f.fig)