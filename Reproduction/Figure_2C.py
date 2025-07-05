#/Users/whyyy/Documents/spatial/protein_correlation_GAT.csv /Users/whyyy/Documents/spatial/protein_correlation_ATT.csv /Users/whyyy/Documents/spatial/protein_correlation_CTP.csv /Users/whyyy/Documents/spatial/protein_correlation_sciPENN.csv /Users/whyyy/Documents/spatial/protein_correlation_Seurat.csv
file_list = [
#'/Users/whyyy/Documents/spatial/Protein_Corr/protein_correlation_FA_DGAT-2.csv',
    './results/Protein_Corr_Leave_One_Out/Spearman_Protein_corr_DGAT.csv',
    './results/Protein_Corr_Leave_One_Out/Spearman_Protein_corr_CTP-net.csv',
    './results/Protein_Corr_Leave_One_Out/Spearman_Protein_corr_scLinear.csv',
    './results/Protein_Corr_Leave_One_Out/Spearman_Protein_corr_sciPENN.csv',
    './results/Protein_Corr_Leave_One_Out/Spearman_Protein_corr_Seurat.csv',
    './results/Protein_Corr_Leave_One_Out/Spearman_mRNA_protein_corr.csv'
]
method_names = {
    'Spearman_Protein_corr_DGAT': 'DGAT',
    'Spearman_Protein_corr_scLinear':'scLinear',
    'Spearman_Protein_corr_CTP-net': 'CTP-net',
     'Spearman_Protein_corr_sciPENN': 'sciPENN',
   'Spearman_Protein_corr_Seurat': 'Seurat v4 (PCA)',
     'Spearman_mRNA_protein_corr': 'mRNA-Protein',
}

def plot_proteinwise_boxplot(file_list, method_names, plot_sample_order=None):
    import os
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load data from CSV files
    data_dict = {}
    for f in file_list:
        key = os.path.basename(f).split('.')[0]
        name = method_names.get(key, key)
        df = pd.read_csv(f, index_col=0)
        data_dict[name] = df

    # Update sample names
    updated_keys = {
        'Tonsil_': 'Tonsil_1',
        'Tonsil_AddOns': 'Tonsil_2',
        'Breast Cancer': 'Breast',
        'MaligMeso_1': 'Meso1',
        'MaligMeso_2': 'Meso_2',
        'Glioblastoma': 'GBM'
    }

    long_list = []
    for method, df in data_dict.items():
        df_melt = (
            df
            .reset_index()
            .melt(id_vars='index', var_name='protein', value_name='corr')
            .rename(columns={'index':'sample'})
        )
        df_melt['sample'] = df_melt['sample'].replace(updated_keys)
        df_melt['method'] = method
        long_list.append(df_melt[['sample','method','corr']])
    df_long = pd.concat(long_list, ignore_index=True)

    samples = df_long['sample'].unique().tolist()
    if plot_sample_order is None:
        plot_sample_order = samples
    methods = list(data_dict.keys())
    puor_colors = ["#8073ac", "#b35806", "#b2abd2", "#e08214", "#fdb863", "#fee0b6"]
    colors = puor_colors[:len(methods)]
    method_colors = {m: colors[i] for i, m in enumerate(methods)}

    width = 0.6
    gap   = 1.0
    positions, data, color_list = [], [], []
    current_pos = 0
    sample_positions = []

    for s in plot_sample_order:
        center = current_pos + (len(methods)-1) * width / 2
        sample_positions.append(center)
        for m in methods:
            arr = df_long.query("sample==@s and method==@m")['corr'].dropna().values
            data.append(arr if arr.size>0 else np.array([np.nan]))
            positions.append(current_pos)
            color_list.append(method_colors[m])
            current_pos += width
        current_pos += gap

    fig, ax = plt.subplots(figsize=(10, 4))
    box = ax.boxplot(data, positions=positions, widths=width, patch_artist=True, showfliers=False)

    for patch, col in zip(box['boxes'], color_list):
        patch.set_facecolor(col)
    for median in box['medians']:
        median.set_color('black')
        median.set_linewidth(1.5)


    ax.set_xticks(sample_positions)
    ax.set_xticklabels(plot_sample_order, rotation=0)

    handles = [plt.Line2D([0],[0], color=method_colors[m], lw=5) for m in methods]
    ax.legend(
        handles, methods,
        title='Method',
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        fontsize=12,
        title_fontsize=12,
        frameon=False
    )
    ax.set_ylabel('Spearman Correlation')
    ax.set_title('Protein-wise Correlation by Sample')
    ax.grid(axis='y', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig('Figure_2C.pdf', dpi=600)

plot_proteinwise_boxplot(file_list, method_names, None)
