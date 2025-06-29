from muon import prot as pt
import scanpy as sc
from collections import defaultdict
import re
import copy
from scipy.sparse import hstack, csr_matrix, issparse
import numpy as np
import pandas as pd
from anndata import AnnData


def qc_control_cytassist(adata, pdata, min_genes=700, max_genes=None, max_mt_pct=None, remove_isotype=True,
                         gene_to_keep_list=None):
    """
    Quality control for CytAssist data.
    :param adata:
    :param pdata:
    :param min_genes:
    :param max_genes:
    :param max_mt_pct:
    :param remove_isotype:
    :param gene_to_keep_list:
    :return:
    """
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    origin_adata = copy.copy(adata)
    # print(adata)
    sc.pp.filter_genes(adata, min_cells=adata.n_obs * 0.025)
    # print(len(adata.var_names))
    gene_list = list(set(gene_to_keep_list) | set(adata.var_names))
    # print(len(gene_list))

    adata = origin_adata[:, gene_list].copy()

    sc.pp.filter_cells(adata, min_genes=min_genes)

    if max_mt_pct is not None:
        adata = adata[adata.obs["pct_counts_mt"] < max_mt_pct]
    if max_genes is not None:
        adata = adata[adata.obs["n_genes_by_counts"] < max_genes]
    if remove_isotype:
        pdata.var["isotype_control"] = (pdata.var_names.str.startswith("mouse_") |
                                        pdata.var_names.str.startswith("rat_") |

                                        pdata.var_names.str.startswith("mouse.") |
                                        pdata.var_names.str.startswith("rat.")
                                        )
        pdata = pdata[:, pdata.var.isotype_control == False]
    pdata = pdata[adata.obs_names, :]
    adata.layers['raw'] = adata.X.copy()
    pdata.layers['raw'] = pdata.X.copy()
    return adata, pdata


def normalize(adata, pdata):
    """
    Normalization for both gene expression and protein data.
    :param adata:
    :param pdata:
    :return:
    """
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    # pdata.X = pdata.layers['raw'].copy()
    pt.pp.clr(pdata)
    # pdata.layers['norm'] = pdata.X.copy()
    return adata, pdata


def clean_protein_names(var_names):
    """
    Clean protein names by replacing hyphens with underscores and parsing the names.
    :param var_names:
    :return:
    """
    # replace - with _
    name_groups = defaultdict(list)
    parsed_names = []

    for vn in var_names:
        vn_clean = vn.replace('-', '_')
        match = re.match(r"^(.*)_(\d+)$", vn_clean)
        if match:
            base, idx = match.groups()
            parsed_names.append((base, int(idx), vn_clean))
            name_groups[base].append(vn_clean)
        else:
            parsed_names.append((vn_clean, None, vn_clean))
            name_groups[vn_clean].append(vn_clean)

    final_names = []
    for base, idx, original in parsed_names:
        if len(name_groups[base]) == 1 and idx is not None:

            final_names.append(base)
        else:
            final_names.append(original)

    return final_names


def QC_train_list(adata_list, pdata_list):
    """
    Quality control for a list of AnnData objects containing gene expression and protein data.
    Finds common genes and proteins across the datasets, performs QC, and normalizes the data.
    :param adata_list:
    :param pdata_list:
    :return:
    """
    gene_list = [set(adata.var_names) for adata in adata_list]
    common_gene = gene_list[0].intersection(*gene_list[1:])
    common_gene = sorted(list(common_gene))
    print(f"Common genes before QC: {len(common_gene)}")

    for pdata in pdata_list:
        pdata.var_names = clean_protein_names(pdata.var_names)
    protein_list = [set(pdata.var_names) for pdata in pdata_list]
    common_protein = protein_list[0].intersection(*protein_list[1:])
    common_protein = sorted(list(common_protein))
    print(f"Common proteins before QC: {len(common_protein)}")

    encoding_gene_list = [
        g for p in common_protein for g in common_gene
        if p.split('_')[0] == g.split('_')[0]
    ]
    print(f"Num of encoding genes: {len(encoding_gene_list)}")

    for adata, pdata, i in zip(adata_list, pdata_list, range(len(adata_list))):
        adata_list[i], pdata_list[i] = qc_control_cytassist(adata_list[i], pdata_list[i], min_genes=700,
                                                            remove_isotype=True, gene_to_keep_list=encoding_gene_list)
        adata_list[i], pdata_list[i] = normalize(adata_list[i], pdata_list[i])

    gene_list = [set(adata.var_names) for adata in adata_list]
    common_gene = gene_list[0].intersection(*gene_list)
    common_gene = sorted(list(common_gene))
    print(f"Common genes after QC: {len(common_gene)}")

    protein_list = [set(pdata.var_names) for pdata in pdata_list]
    common_protein = protein_list[0].intersection(*protein_list)
    common_protein = sorted(list(common_protein))
    print(f"Common proteins after QC: {len(common_protein)}")

    g_filename = f"common_gene_{len(common_gene)}.txt"
    with open(g_filename, "w") as f:
        for gene in common_gene:
            f.write(gene + "\n")
    print(f"Common gene names saved to {g_filename}")

    p_filename = f"common_protein_{len(common_protein)}.txt"
    with open(p_filename, "w") as f:
        for protein in common_protein:
            f.write(protein + "\n")
    print(f"Common protein names saved to {p_filename}")

    return common_gene, common_protein



def fill_genes(test_adata, common_gene):
    """
    Fill missing genes in the test AnnData object with zeros if they are not present in the common gene list.
    :param test_adata:
    :param common_gene:
    :return: adata with filled genes
    """
    existing_genes = set(test_adata.var_names)
    missing_genes = [gene for gene in common_gene if gene not in test_adata.var_names]

    print(f"Sample lost {len(missing_genes)} genes, {len(common_gene)} in total")

    if missing_genes:
        n_obs = test_adata.n_obs
        if issparse(test_adata.X):
            zero_data = csr_matrix((n_obs, len(missing_genes)))
            X_new = hstack([test_adata.X, zero_data]).tocsr()
        else:
            zero_data = np.zeros((n_obs, len(missing_genes)))
            X_new = np.hstack([test_adata.X, zero_data])

        new_var = test_adata.var.copy()

        missing_var = pd.DataFrame(index=missing_genes, columns=new_var.columns)
        new_var = pd.concat([new_var, missing_var], axis=0)
        new_layers = {}
        for k, v in test_adata.layers.items():
            if v is None:
                new_layers[k] = None
            elif issparse(v):

                new_layers[k] = hstack([v, csr_matrix((n_obs, len(missing_genes)))]).tocsr()
            else:
                new_layers[k] = np.hstack([v, zero_data])

        adata = AnnData(
            X=X_new,
            obs=test_adata.obs.copy(),
            var=new_var,
            uns=test_adata.uns.copy(),
            obsm=test_adata.obsm.copy(),
            layers=new_layers,
        )

        test_adata = adata[:, common_gene].copy()
    else:
        test_adata = test_adata[:, common_gene].copy()
    return test_adata


def preprocess_ST(adata):
    """
    Quality control and normalization for spatial transcriptomics data. Due to we selected common genes for the training data, we do not filter genes here.
    :param adata:
    """
    sc.pp.filter_cells(adata, min_genes=700)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)