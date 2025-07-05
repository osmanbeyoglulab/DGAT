from torch_geometric.data import Data, Dataset, HeteroData
from torch_geometric.loader import DataLoader
import pandas as pd

from utils.Graph_utils import MultiGraphDataset, MultiGraphDataset_for_no_protein
from Model.dgat import GATEncoder, Decoder_Protein, Decoder_mRNA
import torch.optim as optim
import torch.nn as nn
from scipy.stats import spearmanr
from tqdm import tqdm
from utils.idk_utils import adata_to_df
from anndata import AnnData
import os

import numpy as np
import torch
import matplotlib.pyplot as plt
import seaborn as sns

alpha, beta, gamma, delta, epsilon, zeta, eta = 5.0, 1.0, 1.0, 3.0, 5.0, 1.0, 1.0
hidden_dim = 1024
dropout_rate = 0.3
batch_size = 1
warmup_epochs = 50
eb_threshold = 0.96


def rmse_loss(pred, target, eps=1e-8):
    return torch.sqrt(nn.MSELoss()(pred, target) + eps)


def train_and_evaluate_fold(train_adata_list, test_adata, train_pdata_list, test_pdata, test_sample_name,
                            processed_data_dir):
    """
    Leave one out testing
    """
    gene_list = [set(adata.var_names) for adata in train_adata_list + [test_adata]]
    common_gene = set.intersection(*gene_list)
    common_gene = sorted(list(common_gene))
    print(f"Common genes: {len(common_gene)}")

    protein_list = [set(pdata.var_names) for pdata in train_pdata_list + [test_pdata]]
    common_protein = set.intersection(*protein_list)
    common_protein = sorted(list(common_protein))
    print(f"Common proteins: {len(common_protein)}")
    with open(f"common_gene_{len(common_gene)}.txt", "w") as f:
        for gene in common_gene:
            f.write(gene + "\n")
    with open(f"common_protein_{len(common_protein)}.txt", "w") as f:
        for protein in common_protein:
            f.write(protein + "\n")

    for i in range(len(train_adata_list)):
        train_adata_list[i] = train_adata_list[i][:, common_gene].copy()
    test_adata = test_adata[:, common_gene].copy()
    for i in range(len(train_pdata_list)):
        train_pdata_list[i] = train_pdata_list[i][:, common_protein].copy()
    test_pdata = test_pdata[:, common_protein].copy()

    dataset = MultiGraphDataset(train_adata_list, train_pdata_list, save_dir=processed_data_dir)
    test_dataset = MultiGraphDataset([test_adata], [test_pdata], save_dir=processed_data_dir)
    train_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=1, shuffle=False)

    if torch.cuda.device_count() >= 2:
        device_mRNA = torch.device("cuda:0")
        device_protein = torch.device("cuda:1")
        print("Using 2 GPUs: mRNA model on cuda:0, protein model on cuda:1")
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        device_mRNA = device
        device_protein = device
        print("Using single device for both models")

    encoder_mRNA = GATEncoder(in_channels=len(common_gene), hidden_dim=hidden_dim, dropout=dropout_rate).to(device_mRNA)
    decoder_mRNA = Decoder_mRNA(hidden_dim, len(common_gene), dropout=0).to(device_mRNA)
    encoder_protein = GATEncoder(in_channels=len(common_protein), hidden_dim=hidden_dim, dropout=dropout_rate).to(
        device_protein)
    decoder_protein = Decoder_Protein(hidden_dim, common_protein, dropout=0).to(device_protein)

    optimizer_enc_mRNA = optim.Adam(encoder_mRNA.parameters(), lr=5e-4, weight_decay=2e-5)
    optimizer_dec_mRNA = optim.Adam(decoder_mRNA.parameters(), lr=1e-4, weight_decay=1e-5)
    optimizer_enc_protein = optim.Adam(encoder_protein.parameters(), lr=5e-4, weight_decay=2e-5)
    optimizer_dec_protein = optim.Adam(decoder_protein.parameters(), lr=1e-4, weight_decay=1e-5)
    mse_loss_fn = nn.MSELoss()

    epochs = 100
    eval_interval = 1
    train_losses = []
    test_losses = []
    test_cors = []


    D = sum(p.numel() for p in encoder_mRNA.parameters()) + \
        sum(p.numel() for p in decoder_mRNA.parameters()) + \
        sum(p.numel() for p in encoder_protein.parameters()) + \
        sum(p.numel() for p in decoder_protein.parameters())

    evidence_counters = []

    for epoch in range(1, epochs + 1):
        encoder_mRNA.train()
        decoder_mRNA.train()
        encoder_protein.train()
        decoder_protein.train()
        total_loss = 0.0
        valid_batches = 0

        for data in tqdm(train_loader, desc=f"Epoch {epoch}/{epochs}", leave=False):
            mRNA_batch = data['mRNA'].x.to(device_mRNA)
            mRNA_edge_index = data[('mRNA', 'mRNA_knn', 'mRNA')].edge_index.to(device_mRNA)
            mRNA_edge_attr = data[('mRNA', 'mRNA_knn', 'mRNA')].edge_attr.to(device_mRNA)
            protein_batch = data['protein'].x.to(device_protein)
            protein_edge_index = data[('protein', 'protein_knn', 'protein')].edge_index.to(device_protein)
            protein_edge_attr = data[('protein', 'protein_knn', 'protein')].edge_attr.to(device_protein)

            optimizer_enc_mRNA.zero_grad()
            optimizer_dec_mRNA.zero_grad()
            optimizer_enc_protein.zero_grad()
            optimizer_dec_protein.zero_grad()

            try:
                z_mRNA = encoder_mRNA(mRNA_batch, mRNA_edge_index, mRNA_edge_attr)
                z_protein = encoder_protein(protein_batch, protein_edge_index, protein_edge_attr)
            except Exception as e:
                print(f"Error in encoder forward pass: {e}")
                continue

            mRNA_recon = decoder_mRNA(z_mRNA)
            protein_recon_from_protein = decoder_protein(z_protein)

            with torch.no_grad():
                z_protein_on_mRNA = z_protein.to(device_mRNA, non_blocking=True)

            loss_mRNA_recon = rmse_loss(mRNA_recon, mRNA_batch)
            loss_protein_recon = rmse_loss(protein_recon_from_protein, protein_batch)
            loss_align = mse_loss_fn(z_mRNA, z_protein_on_mRNA)
            loss_protein_pred = rmse_loss(decoder_protein(z_mRNA.to(device_protein, non_blocking=True)), protein_batch)
            loss_mRNA_pred = rmse_loss(decoder_mRNA(z_protein_on_mRNA), mRNA_batch)

            loss_protein_recon = loss_protein_recon.to(device_mRNA, non_blocking=True)
            loss_protein_pred = loss_protein_pred.to(device_mRNA, non_blocking=True)

            if loss_align < 0.015:
                effective_gamma = 0.0
            else:
                effective_gamma = gamma

            if loss_protein_pred < 0.015:
                effective_delta = 0.0
            else:
                effective_delta = delta

            if loss_protein_recon < 0.015:
                effective_beta = 0.0
            else:
                effective_beta = beta

            total_loss_batch = (alpha * loss_mRNA_recon + effective_beta * loss_protein_recon +
                                effective_gamma * loss_align + effective_delta * loss_protein_pred + eta * loss_mRNA_pred)

            total_loss_batch.backward()
            torch.nn.utils.clip_grad_norm_(encoder_mRNA.parameters(), max_norm=1.0)
            torch.nn.utils.clip_grad_norm_(decoder_mRNA.parameters(), max_norm=1.0)
            torch.nn.utils.clip_grad_norm_(encoder_protein.parameters(), max_norm=1.0)
            torch.nn.utils.clip_grad_norm_(decoder_protein.parameters(), max_norm=1.0)
            optimizer_enc_mRNA.step()
            optimizer_dec_mRNA.step()
            optimizer_enc_protein.step()
            optimizer_dec_protein.step()

            total_loss += total_loss_batch.item()
            valid_batches += 1

        avg_loss = total_loss / valid_batches if valid_batches > 0 else 0
        print(f"Epoch [{epoch}/{epochs}] Total Loss: {avg_loss:.4f} | "
              f"mRNA Recon: {loss_mRNA_recon:.4f} | "
              f"Protein Recon: {loss_protein_recon:.4f} | "
              f"Alignment: {loss_align:.4f} | "
              f"Protein Pred: {loss_protein_pred:.4f} | "
              f"mRNA Pred: {loss_mRNA_pred:.4f}")
        train_losses.append(avg_loss)
        torch.cuda.empty_cache()
        if epoch % eval_interval == 0:
            test_loss, test_cor = test_data_spearman_only(train_loader, test_loader, encoder_mRNA, decoder_mRNA,
                                                          encoder_protein, decoder_protein, device_mRNA, device_protein,
                                                          mse_loss_fn, test_adata, test_pdata, common_protein,
                                                          test_sample_name)
            test_losses.append(test_loss)
            test_cors.append(test_cor)
        torch.cuda.empty_cache()

        # lr decay
        if epoch % 10 == 0:
            for optimizer in [optimizer_enc_mRNA, optimizer_dec_mRNA, optimizer_enc_protein, optimizer_dec_protein]:
                for param_group in optimizer.param_groups:
                    param_group['lr'] *= 0.8
            print(f"Epoch {epoch}: Learning rate reduced.")

        # ——— EB‐criterion computation ———
        sum_term = 0.0
        for opt in (optimizer_enc_mRNA, optimizer_dec_mRNA, optimizer_enc_protein, optimizer_dec_protein):
            for grp in opt.param_groups:
                for p in grp['params']:
                    if p.grad is None:
                        continue
                    st = opt.state[p]
                    m_t = st.get('exp_avg')
                    v_t = st.get('exp_avg_sq')
                    if m_t is None or v_t is None:
                        continue
                    g = p.grad.data
                    var = v_t - m_t.pow(2) + 1e-8
                    sum_term += (g.pow(2) / var).sum().item()
        evidence = 1.0 - (batch_size / D) * sum_term

        print(f"Epoch {epoch}: EB evidence = {evidence:.4f}")
        # maintain sliding window of last 10 evidences
        evidence_counters.append(evidence)
        if len(evidence_counters) > 10:
            evidence_counters.pop(0)
        mean_evidence = sum(evidence_counters) / len(evidence_counters)

        # early stop if after warmup and mean_evidence exceeds threshold
        if epoch > warmup_epochs and mean_evidence > eb_threshold:
            print(f"--> EB early stopping at epoch {epoch} "
                  f"(mean_evidence={mean_evidence:.4f}, threshold={eb_threshold})")
            break

    final_test_loss, final_cor = test_data_spearman_only(train_loader, test_loader, encoder_mRNA, decoder_mRNA,
                                                         encoder_protein, decoder_protein, device_mRNA, device_protein,
                                                         mse_loss_fn, test_adata, test_pdata, common_protein,
                                                         test_sample_name, if_tt=True)
    model_components = {
        'mRNA_encoder': encoder_mRNA,
        'mRNA_decoder': decoder_mRNA,
        'protein_encoder': encoder_protein,
        'protein_decoder': decoder_protein

    }

    return model_components


def val_cor(adata_1, pdata_1_df, protein_names=None, if_tt=False, title=''):
    global common_protein
    # get_protein_expressions
    a_1_p = get_activity(adata_1, key='protein_predict', protein_names=protein_names)
    df_a_1_p = adata_to_df(a_1_p)
    df_pdata_1 = pdata_1_df.copy()
    common_obs = list(set(a_1_p.obs_names) & set(df_pdata_1.index))
    a_1_p = a_1_p[common_obs, :].copy()
    df_a_1_p = df_a_1_p.loc[common_obs, :]
    df_pdata_1 = df_pdata_1.loc[common_obs, :]

    common_protein_chk = list(set(df_a_1_p.columns) & set(df_pdata_1.columns))
    df_a_1_p = df_a_1_p[common_protein_chk]
    df_pdata_1 = df_pdata_1[common_protein_chk]

    corrs = df_a_1_p.corrwith(df_pdata_1, method='pearson')
    corrs = corrs.sort_values(ascending=False)

    row_corrs = df_a_1_p.T.corrwith(df_pdata_1.T, method='pearson')
    row_corrs = df_a_1_p.apply(lambda row: row.corr(df_pdata_1.loc[row.name], method='pearson'), axis=1)
    os.makedirs('./Corrs_All_Spots_My_model', exist_ok=True)
    row_corrs.to_csv(f'./Corrs_All_Spots_My_model/Corrs_All_Spots_{title}_My_model.csv', index=True)

    print('Average Cell Pearson Cor in test data: ', row_corrs.mean())
    row_corrs = df_a_1_p.T.corrwith(df_pdata_1.T, method='spearman')
    row_corrs = df_a_1_p.apply(lambda row: row.corr(df_pdata_1.loc[row.name], method='spearman'), axis=1)
    os.makedirs('./Corrs_All_Spots_My_model', exist_ok=True)
    row_corrs.to_csv(f'./Corrs_All_Spots_My_model/Corrs_All_Spots_{title}_My_model.csv', index=True)

    print('Average Cell Spearman Cor in test data: ', row_corrs.mean())
    mean_corr = corrs.mean()
    print(f"Average Protein Cor in test data: {mean_corr:.4f}\n")
    rmse = np.sqrt(np.mean((df_a_1_p - df_pdata_1).values ** 2))
    print(f"Average RMSE in test data: {rmse:.4f}\n")

    if not if_tt:
        return

    order_list = protein_names
    new = pd.DataFrame(corrs, columns=['correlation']).reset_index()
    new.columns = ['protein_names', 'correlation']
    new = new.set_index('protein_names')
    new = new.loc[order_list]

    if not title:
        title = 'sample'
    df_save = new.copy()
    df_save.rename(columns={'correlation': title}, inplace=True)
    df_save = df_save.reset_index()

    plt.figure(figsize=(15, 1))
    sns.heatmap(new.T, cmap='coolwarm', annot=True, fmt=".2f", cbar=False, linewidths=5)
    plt.title(f'Correlation between predicted protein expression and protein expression {title}')
    # plt.savefig(f'./GCN_result/Test_Cor_{title}_My_model.png')
    plt.show()

    num_proteins = len(order_list)
    n_cols = 5
    n_rows = int(np.ceil(num_proteins / n_cols))
    fig = plt.figure(figsize=(n_cols * 4, n_rows * 4))
    fig.suptitle(f'Scatter Plot of Predicted vs True Expression for {title}', fontsize=16)
    for idx, protein in enumerate(order_list):
        if protein not in df_a_1_p.columns:
            continue
        pred_values = df_a_1_p[protein]
        true_values = df_pdata_1[protein]
        min_val = min(pred_values.min(), true_values.min())
        max_val = max(pred_values.max(), true_values.max())
        margin = (max_val - min_val) * 0.05
        axis_min = min_val - margin
        axis_max = max_val + margin
        corr_val = pred_values.corr(true_values)
        plt.subplot(n_rows, n_cols, idx + 1)
        sns.scatterplot(x=pred_values, y=true_values, s=20, color='blue')
        slope, intercept = np.polyfit(pred_values, true_values, 1)
        x_fit = np.linspace(axis_min, axis_max, 100)
        y_fit = slope * x_fit + intercept
        plt.plot(x_fit, y_fit, 'r--', lw=1)
        plt.xlim(axis_min, axis_max)
        plt.ylim(axis_min, axis_max)
        plt.xlabel('Predicted')
        plt.ylabel('True')
        plt.title(protein)
        plt.text(0.05, 0.95, f'r={corr_val:.2f}', transform=plt.gca().transAxes,
                 fontsize=10, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    # plt.savefig(f'./Scatter_{title}_my_Model.png')
    plt.show()


def train(train_adata_list, train_pdata_list, processed_data_dir):
    gene_list = [set(adata.var_names) for adata in train_adata_list]
    common_gene = set.intersection(*gene_list)
    common_gene = sorted(list(common_gene))
    print(f"Common genes: {len(common_gene)}")

    protein_list = [set(pdata.var_names) for pdata in train_pdata_list]
    common_protein = set.intersection(*protein_list)
    common_protein = sorted(list(common_protein))
    print(f"Common proteins: {len(common_protein)}")
    with open(f"common_gene_{len(common_gene)}.txt", "w") as f:
        for gene in common_gene:
            f.write(gene + "\n")
    with open(f"common_protein_{len(common_protein)}.txt", "w") as f:
        for protein in common_protein:
            f.write(protein + "\n")

    for i in range(len(train_adata_list)):
        train_adata_list[i] = train_adata_list[i][:, common_gene].copy()
    # test_adata = test_adata[:, common_gene].copy()
    for i in range(len(train_pdata_list)):
        train_pdata_list[i] = train_pdata_list[i][:, common_protein].copy()
    # test_pdata = test_pdata[:, common_protein].copy()

    dataset = MultiGraphDataset(train_adata_list, train_pdata_list, save_dir=processed_data_dir)
    # test_dataset = MultiGraphDataset([test_adata], [test_pdata], save_dir=processed_data_dir)
    train_loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    # test_loader = DataLoader(test_dataset, batch_size=1, shuffle=False)

    if torch.cuda.device_count() >= 2:
        device_mRNA = torch.device("cuda:0")
        device_protein = torch.device("cuda:1")
        print("Using 2 GPUs: mRNA model on cuda:0, protein model on cuda:1")
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        device_mRNA = device
        device_protein = device
        print("Using single device for both models")

    encoder_mRNA = GATEncoder(in_channels=len(common_gene), hidden_dim=hidden_dim, dropout=dropout_rate).to(device_mRNA)
    decoder_mRNA = Decoder_mRNA(hidden_dim, len(common_gene), dropout=0).to(device_mRNA)
    encoder_protein = GATEncoder(in_channels=len(common_protein), hidden_dim=hidden_dim, dropout=dropout_rate).to(
        device_protein)
    decoder_protein = Decoder_Protein(hidden_dim, common_protein, dropout=0).to(device_protein)

    optimizer_enc_mRNA = optim.Adam(encoder_mRNA.parameters(), lr=5e-4, weight_decay=2e-5)
    optimizer_dec_mRNA = optim.Adam(decoder_mRNA.parameters(), lr=1e-4, weight_decay=1e-5)
    optimizer_enc_protein = optim.Adam(encoder_protein.parameters(), lr=5e-4, weight_decay=2e-5)
    optimizer_dec_protein = optim.Adam(decoder_protein.parameters(), lr=1e-4, weight_decay=1e-5)
    mse_loss_fn = nn.MSELoss()

    epochs = 100
    eval_interval = 1
    train_losses = []



    D = sum(p.numel() for p in encoder_mRNA.parameters()) + \
        sum(p.numel() for p in decoder_mRNA.parameters()) + \
        sum(p.numel() for p in encoder_protein.parameters()) + \
        sum(p.numel() for p in decoder_protein.parameters())

    evidence_counters = []

    for epoch in range(1, epochs + 1):
        encoder_mRNA.train()
        decoder_mRNA.train()
        encoder_protein.train()
        decoder_protein.train()
        total_loss = 0.0
        valid_batches = 0

        for data in tqdm(train_loader, desc=f"Epoch {epoch}/{epochs}", leave=False):
            mRNA_batch = data['mRNA'].x.to(device_mRNA)
            mRNA_edge_index = data[('mRNA', 'mRNA_knn', 'mRNA')].edge_index.to(device_mRNA)
            mRNA_edge_attr = data[('mRNA', 'mRNA_knn', 'mRNA')].edge_attr.to(device_mRNA)
            protein_batch = data['protein'].x.to(device_protein)
            protein_edge_index = data[('protein', 'protein_knn', 'protein')].edge_index.to(device_protein)
            protein_edge_attr = data[('protein', 'protein_knn', 'protein')].edge_attr.to(device_protein)

            optimizer_enc_mRNA.zero_grad()
            optimizer_dec_mRNA.zero_grad()
            optimizer_enc_protein.zero_grad()
            optimizer_dec_protein.zero_grad()

            try:
                z_mRNA = encoder_mRNA(mRNA_batch, mRNA_edge_index, mRNA_edge_attr)
                z_protein = encoder_protein(protein_batch, protein_edge_index, protein_edge_attr)
            except Exception as e:
                print(f"Error in encoder forward pass: {e}")
                continue

            mRNA_recon = decoder_mRNA(z_mRNA)
            protein_recon_from_protein = decoder_protein(z_protein)

            with torch.no_grad():
                z_protein_on_mRNA = z_protein.to(device_mRNA, non_blocking=True)

            loss_mRNA_recon = rmse_loss(mRNA_recon, mRNA_batch)
            loss_protein_recon = rmse_loss(protein_recon_from_protein, protein_batch)
            loss_align = mse_loss_fn(z_mRNA, z_protein_on_mRNA)
            loss_protein_pred = rmse_loss(decoder_protein(z_mRNA.to(device_protein, non_blocking=True)), protein_batch)
            loss_mRNA_pred = rmse_loss(decoder_mRNA(z_protein_on_mRNA), mRNA_batch)

            loss_protein_recon = loss_protein_recon.to(device_mRNA, non_blocking=True)
            loss_protein_pred = loss_protein_pred.to(device_mRNA, non_blocking=True)

            if loss_align < 0.015:
                effective_gamma = 0.0
            else:
                effective_gamma = gamma

            if loss_protein_pred < 0.015:
                effective_delta = 0.0
            else:
                effective_delta = delta

            if loss_protein_recon < 0.015:
                effective_beta = 0.0
            else:
                effective_beta = beta

            total_loss_batch = (alpha * loss_mRNA_recon + effective_beta * loss_protein_recon +
                                effective_gamma * loss_align + effective_delta * loss_protein_pred + eta * loss_mRNA_pred)

            total_loss_batch.backward()
            torch.nn.utils.clip_grad_norm_(encoder_mRNA.parameters(), max_norm=1.0)
            torch.nn.utils.clip_grad_norm_(decoder_mRNA.parameters(), max_norm=1.0)
            torch.nn.utils.clip_grad_norm_(encoder_protein.parameters(), max_norm=1.0)
            torch.nn.utils.clip_grad_norm_(decoder_protein.parameters(), max_norm=1.0)
            optimizer_enc_mRNA.step()
            optimizer_dec_mRNA.step()
            optimizer_enc_protein.step()
            optimizer_dec_protein.step()

            total_loss += total_loss_batch.item()
            valid_batches += 1

        avg_loss = total_loss / valid_batches if valid_batches > 0 else 0
        print(f"Epoch [{epoch}/{epochs}] Total Loss: {avg_loss:.4f} | "
              f"mRNA Recon: {loss_mRNA_recon:.4f} | "
              f"Protein Recon: {loss_protein_recon:.4f} | "
              f"Alignment: {loss_align:.4f} | "
              f"Protein Pred: {loss_protein_pred:.4f} | "
              f"mRNA Pred: {loss_mRNA_pred:.4f}")
        train_losses.append(avg_loss)
        torch.cuda.empty_cache()
        torch.cuda.empty_cache()

        # lr decay
        if epoch % 10 == 0:
            for optimizer in [optimizer_enc_mRNA, optimizer_dec_mRNA, optimizer_enc_protein, optimizer_dec_protein]:
                for param_group in optimizer.param_groups:
                    param_group['lr'] *= 0.8
            print(f"Epoch {epoch}: Learning rate reduced.")

        # ——— EB‐criterion computation ———
        sum_term = 0.0
        for opt in (optimizer_enc_mRNA, optimizer_dec_mRNA, optimizer_enc_protein, optimizer_dec_protein):
            for grp in opt.param_groups:
                for p in grp['params']:
                    if p.grad is None:
                        continue
                    st = opt.state[p]
                    m_t = st.get('exp_avg')
                    v_t = st.get('exp_avg_sq')
                    if m_t is None or v_t is None:
                        continue
                    g = p.grad.data
                    var = v_t - m_t.pow(2) + 1e-8
                    sum_term += (g.pow(2) / var).sum().item()
        evidence = 1.0 - (batch_size / D) * sum_term

        print(f"Epoch {epoch}: EB evidence = {evidence:.4f}")
        # maintain sliding window of last N evidences
        evidence_counters.append(evidence)
        if len(evidence_counters) > 10:
            evidence_counters.pop(0)
        mean_evidence = sum(evidence_counters) / len(evidence_counters)

        # early stop if after warmup and mean_evidence exceeds threshold
        if epoch > warmup_epochs and mean_evidence > eb_threshold:
            print(f"--> EB early stopping at epoch {epoch} "
                  f"(mean_evidence={mean_evidence:.4f}, threshold={eb_threshold})")
            break
    model_components = {
        'mRNA_encoder': encoder_mRNA,
        'mRNA_decoder': decoder_mRNA,
        'protein_encoder': encoder_protein,
        'protein_decoder': decoder_protein

    }

    return model_components


def test_data_spearman_only(
        train_loader, test_loader,
        encoder_mRNA, decoder_mRNA,
        encoder_protein, decoder_protein,
        device_mRNA, device_protein,
        mse_loss_fn, test_adata, test_pdata,
        common_protein, test_sample_name,
        if_tt=False
):
    encoder_mRNA.eval()
    decoder_mRNA.eval()
    encoder_protein.eval()
    decoder_protein.eval()

    total_loss = 0.0
    valid_batches = 0
    all_protein_pred = []
    all_protein_true = []

    with torch.no_grad():
        for data in test_loader:
            # mRNA graph
            data = data.to(device_mRNA)
            mRNA_test = data['mRNA'].x.to(device_mRNA)
            mRNA_edge_index = data[('mRNA', 'mRNA_knn', 'mRNA')].edge_index.to(device_mRNA)
            mRNA_edge_attr = data[('mRNA', 'mRNA_knn', 'mRNA')].edge_attr.to(device_mRNA)
            z_mRNA = encoder_mRNA(mRNA_test, mRNA_edge_index, mRNA_edge_attr)
            mRNA_recon = decoder_mRNA(z_mRNA)

            # protein graph
            protein_test = data['protein'].x.to(device_protein)
            protein_edge_index = data[('protein', 'protein_knn', 'protein')].edge_index.to(device_protein)
            protein_edge_attr = data[('protein', 'protein_knn', 'protein')].edge_attr.to(device_protein)
            z_protein = encoder_protein(protein_test, protein_edge_index, protein_edge_attr)
            protein_recon = decoder_protein(z_protein)

            # cross-modal prediction
            protein_recon_from_mRNA = decoder_protein(z_mRNA.to(device_protein))
            mRNA_recon_from_protein = decoder_mRNA(z_protein.to(device_mRNA))

            # losses
            loss_mRNA_recon = rmse_loss(mRNA_recon, mRNA_test)
            loss_protein_recon = rmse_loss(protein_recon, protein_test)
            loss_align = mse_loss_fn(z_mRNA, z_protein.to(device_mRNA))
            loss_protein_pred = rmse_loss(protein_recon_from_mRNA, protein_test)
            loss_mRNA_pred = rmse_loss(mRNA_recon_from_protein, mRNA_test)

            # total loss
            loss_protein_recon = loss_protein_recon.to(device_mRNA)
            loss_protein_pred = loss_protein_pred.to(device_mRNA)
            batch_loss = (alpha * loss_mRNA_recon + beta * loss_protein_recon +
                          gamma * loss_align + delta * loss_protein_pred +
                          epsilon * loss_mRNA_pred)
            total_loss += batch_loss.item()
            valid_batches += 1

            all_protein_pred.append(protein_recon_from_mRNA.cpu().numpy())
            all_protein_true.append(protein_test.cpu().numpy())

    avg_loss = total_loss / valid_batches if valid_batches > 0 else 0
    all_protein_pred = np.concatenate(all_protein_pred, axis=0)
    all_protein_true = np.concatenate(all_protein_true, axis=0)

    spearman_corrs_pro = [
        spearmanr(all_protein_pred[:, i], all_protein_true[:, i])[0]
        for i in range(all_protein_pred.shape[1])
    ]
    avg_spearman_Pro = np.nanmean(spearman_corrs_pro)

    spearman_corrs_cells = []
    low_var_count = 0
    for i in range(all_protein_pred.shape[0]):
        pred_vec = all_protein_pred[i, :]
        true_vec = all_protein_true[i, :]
        if np.std(pred_vec) < 1e-6 or np.std(true_vec) < 1e-6:
            low_var_count += 1
            continue
        sp_corr, _ = spearmanr(pred_vec, true_vec)
        spearman_corrs_cells.append(sp_corr)

    spearman_corrs_cells = np.array(spearman_corrs_cells)
    avg_spearman_cells = np.nanmean(spearman_corrs_cells)

    print(f"Test Loss: {avg_loss:.4f}")
    print("Protein-wise Spearman: Mean = {:.4f}, Std = {:.4f}, Median = {:.4f}".format(
        np.nanmean(spearman_corrs_pro),
        np.nanstd(spearman_corrs_pro),
        np.nanmedian(spearman_corrs_pro)
    ))
    print("Cell-wise Spearman: Mean = {:.4f}, Std = {:.4f}, Median = {:.4f}".format(
        np.nanmean(spearman_corrs_cells),
        np.nanstd(spearman_corrs_cells),
        np.nanmedian(spearman_corrs_cells)
    ))
    print("Number of cells with low variation:", low_var_count)

    if if_tt:
        test_adata.obsm['protein_predict'] = all_protein_pred
        test_pdata_df = adata_to_df(test_pdata)
        val_cor(test_adata, test_pdata_df, protein_names=list(common_protein), if_tt=if_tt, title=test_sample_name)

    return avg_loss, avg_spearman_Pro


def get_activity(adata, key='', protein_names=None):
    s = key[:2]
    adata_pt = AnnData(
        X=adata.obsm[key],
        obs=adata.obs,
        obsm={name: obj for (name, obj) in adata.obsm.items() if s not in name},
        layers={name: obj for (name, obj) in adata.obsm.items() if s in name},
    )
    if protein_names is not None:
        adata_pt.var_names = protein_names
    adata_pt.uns = adata.uns
    return adata_pt


def predict_data(ae_mRNA, decoder_protein, test_adata, predict_loader, device):
    ae_mRNA.eval()
    decoder_protein.eval()
    with torch.no_grad():
        for data in predict_loader:
            data = data.to(device)
            mRNA_test = data['mRNA'].x.to(device)
            try:
                z_mRNA_test = ae_mRNA(mRNA_test,
                                      data[('mRNA', 'mRNA_knn', 'mRNA')].edge_index.to(device),
                                      data[('mRNA', 'mRNA_knn', 'mRNA')].edge_attr.to(device))
            except Exception as e:
                print(f"Error in test mRNA forward pass: {e}")
                continue
            # device = next(decoder_protein.parameters()).device
            # print(f"Model is on: {device}")
            z_mRNA_test = z_mRNA_test.to(device)
            decoder_protein.to(device)
            # decoder_mRNA.to(device_mRNA)
            # mRNA_recon = decoder_mRNA(z_mRNA_test)
            protein_pred = decoder_protein(z_mRNA_test)  # .to(device_protein))
            protein_pred_all = protein_pred.cpu().numpy()
            # mRNA_recon_all = mRNA_recon.cpu().numpy()
            test_adata.obsm['protein_predict'] = protein_pred_all
            # test_adata.obsm['mRNA_recon'] = mRNA_recon_all


def protein_predict(adata, common_gene, common_protein, model_save_dir, pyg_data_dir):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    dir_path = f"{model_save_dir}/{len(common_gene)}_gene_{len(common_protein)}_protein"
    encoder_mRNA = GATEncoder(in_channels=len(common_gene), hidden_dim=hidden_dim, dropout=dropout_rate)
    decoder_protein = Decoder_Protein(hidden_dim, common_protein)

    encoder_mRNA.load_state_dict(
        torch.load(os.path.join(dir_path, f'encoder_mRNA.pth'), map_location=torch.device('cpu')))
    decoder_protein.load_state_dict(
        torch.load(os.path.join(dir_path, f'decoder_protein.pth'), map_location=torch.device('cpu')))

    encoder_mRNA.eval()
    decoder_protein.eval()
    encoder_mRNA.to(device)
    decoder_protein.to(device)
    test_dataset = MultiGraphDataset_for_no_protein([adata], device=device, save_dir=pyg_data_dir)
    predict_loader = DataLoader(test_dataset, batch_size=1, shuffle=False)
    predict_data(encoder_mRNA, decoder_protein, adata, predict_loader, device)

    return get_activity(adata, key='protein_predict', protein_names=common_protein)
