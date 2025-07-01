from torch_geometric.data import Data, Dataset, HeteroData
import os
import torch
from sklearn.neighbors import NearestNeighbors
from utils.idk_utils import adata_to_df
import numpy as np
from sklearn.decomposition import PCA


def build_knn_adj(features, k, apply_pca=False, variance=0.85):
    """
    Build a k-nearest neighbors adjacency matrix from the given features.
    """
    if apply_pca and features.shape[1] > 1500:
        pca = PCA(n_components=variance, svd_solver='full')
        features = pca.fit_transform(features)
    nbrs = NearestNeighbors(n_neighbors=k, algorithm='ball_tree').fit(features)
    distances, indices = nbrs.kneighbors(features)
    num_nodes = features.shape[0]
    adjacency = torch.zeros((num_nodes, num_nodes), dtype=torch.float32)
    for i, neighbors in enumerate(indices):
        for neighbor in neighbors:
            adjacency[i, neighbor] = 1.0
    adjacency += torch.eye(num_nodes) * 1
    return adjacency


def preprocess_data(adata, pdata_df):
    X_mRNA = adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X
    X_protein = pdata_df.values
    spatial = adata.obsm['spatial']
    return X_mRNA, X_protein, spatial


def adjacency_to_edge_index(adjacency):
    """
    Convert an adjacency matrix to edge index and edge weight tensors.
    """
    edge_index = adjacency.nonzero(as_tuple=False).t()
    edge_weight = adjacency[edge_index[0], edge_index[1]]
    return edge_index.long(), edge_weight


def create_pyg_data(X_mRNA, X_protein, spatial, device='cpu'):
    """
    Create a PyTorch Geometric data object from mRNA and protein expression data, along with spatial coordinates.
    """
    adjacency_spatial = build_knn_adj(spatial, k=6, apply_pca=False)
    adjacency_mRNA = build_knn_adj(X_mRNA, k=10, apply_pca=True, variance=0.85)
    adjacency_protein = build_knn_adj(X_protein, k=10, apply_pca=False)
    adj_mRNA = adjacency_spatial + adjacency_mRNA
    adj_protein = adjacency_spatial + adjacency_protein
    edge_index_mRNA, edge_weight_mRNA = adjacency_to_edge_index(adj_mRNA)
    edge_index_protein, edge_weight_protein = adjacency_to_edge_index(adj_protein)
    if edge_index_mRNA.numel() == 0:
        raise ValueError("edge_index_mRNA is empty.")
    if edge_index_protein.numel() == 0:
        raise ValueError("edge_index_protein is empty.")
    data = HeteroData()
    data['mRNA'].x = torch.tensor(X_mRNA, dtype=torch.float32)
    data['protein'].x = torch.tensor(X_protein, dtype=torch.float32)
    data['mRNA'].pos = torch.tensor(spatial, dtype=torch.float32)
    data['protein'].pos = torch.tensor(spatial, dtype=torch.float32)
    data['mRNA', 'mRNA_knn', 'mRNA'].edge_index = edge_index_mRNA
    data['mRNA', 'mRNA_knn', 'mRNA'].edge_attr = edge_weight_mRNA
    data['protein', 'protein_knn', 'protein'].edge_index = edge_index_protein
    data['protein', 'protein_knn', 'protein'].edge_attr = edge_weight_protein
    return data


def save_pyg_data(pyg_data, filepath):
    torch.save(pyg_data, filepath)


def load_pyg_data(filepath):
    return torch.load(filepath)


class MultiGraphDataset(Dataset):
    def __init__(self, adata_list, pdata_list, device='cuda', save_dir='processed_data'):
        super(MultiGraphDataset, self).__init__()
        self.data_list = []
        self.save_dir = save_dir
        os.makedirs(self.save_dir, exist_ok=True)
        print("Creating or loading dataset")
        for idx, (adata, pdata) in enumerate(zip(adata_list, pdata_list)):
            try:
                sample_name = list(adata.uns['spatial'].keys())[0]
            except:
                sample_name = adata.uns['name']

            filepath = os.path.join(self.save_dir,
                                    f"{sample_name}_{adata.shape[1]}_new_spatial=6_1hop_new_filtering.pth")
            if os.path.exists(filepath):
                print(f"Loading preprocessed data for sample '{sample_name}' from '{filepath}'")
                pyg_data = load_pyg_data(filepath)
            else:
                print(f"Creating and saving preprocessed data for sample '{sample_name}'")
                pdata_df = adata_to_df(pdata)
                X_mRNA, X_protein, spatial = preprocess_data(adata, pdata_df)
                pyg_data = create_pyg_data(X_mRNA, X_protein, spatial, device=device)
                save_pyg_data(pyg_data, filepath)
            self.data_list.append(pyg_data)
        print("Dataset ready")

    def len(self):
        return len(self.data_list)

    def __len__(self):
        return len(self.data_list)

    def get(self, idx):
        return self.data_list[idx]


# test_pdata = adata.copy()
def preprocess_data_no_protein(adata_p):
    X_mRNA = adata_p.X.toarray() if not isinstance(adata_p.X, np.ndarray) else adata_p.X
    # X_protein = pdata_df.values
    spatial = adata_p.obsm['spatial']
    return X_mRNA, spatial


def create_pyg_data_no_protein(X_mRNA, spatial, device='cpu'):
    adjacency_spatial = build_knn_adj(spatial, k=6, apply_pca=False, variance=0.85)
    adjacency_mRNA = build_knn_adj(X_mRNA, k=10, apply_pca=True, variance=0.85)
    # adjacency_protein = build_knn_adj(X_protein, lof_proteins, k=10, apply_pca=False, variance=0.85)

    adj_mRNA = adjacency_spatial + adjacency_mRNA
    # adj_protein = adjacency_spatial + adjacency_protein

    edge_index_mRNA, edge_weight_mRNA = adjacency_to_edge_index(adj_mRNA)
    # edge_index_protein, edge_weight_protein = adjacency_to_edge_index(adj_protein)

    if edge_index_mRNA.numel() == 0:
        raise ValueError("edge_index_mRNA is empty.")
    # if edge_index_protein.numel() == 0:
    #     raise ValueError("edge_index_protein is empty.")

    data = HeteroData()
    data['mRNA'].x = torch.tensor(X_mRNA, dtype=torch.float32)
    # data['protein'].x = torch.tensor(X_protein, dtype=torch.float32)

    data['mRNA'].pos = torch.tensor(spatial, dtype=torch.float32)
    # data['protein'].pos = torch.tensor(spatial, dtype=torch.float32)

    data['mRNA', 'mRNA_knn', 'mRNA'].edge_index = edge_index_mRNA
    data['mRNA', 'mRNA_knn', 'mRNA'].edge_attr = edge_weight_mRNA
    # data['protein', 'protein_knn', 'protein'].edge_index = edge_index_protein
    # data['protein', 'protein_knn', 'protein'].edge_attr = edge_weight_protein

    return data


class MultiGraphDataset_for_no_protein(Dataset):
    def __init__(self, adata_list, device='cuda', save_dir='processed_data'):
        super(MultiGraphDataset_for_no_protein, self).__init__()
        self.data_list = []
        self.save_dir = save_dir
        os.makedirs(self.save_dir, exist_ok=True)
        print("Creating or loading dataset")

        for idx, (adata_p) in enumerate(adata_list):
            try:
                sample_name = list(adata_p.uns['spatial'].keys())[0]
            except:
                sample_name = adata_p.uns['name']
            filepath = os.path.join(self.save_dir, f"{sample_name}_{adata_p.shape[1]}_spatial=6_testing.pth")

            if os.path.exists(filepath):
                print(f"Loading preprocessed data for sample '{sample_name}' from '{filepath}'")
                pyg_data = load_pyg_data(filepath)
            else:
                print(f"Creating and saving preprocessed data for sample '{sample_name}'")
                X_mRNA, spatial = preprocess_data_no_protein(adata_p)

                pyg_data = create_pyg_data_no_protein(X_mRNA, spatial, device=device)
                save_pyg_data(pyg_data, filepath)

            self.data_list.append(pyg_data)
        print("Dataset ready")

    def len(self):
        return len(self.data_list)

    def __len__(self):
        return len(self.data_list)

    def get(self, idx):
        return self.data_list[idx]


processed_data_dir = '/ix/hosmanbeyoglu/haw309/new_pyg_data'



