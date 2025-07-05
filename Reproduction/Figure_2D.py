import pandas as pd
import matplotlib.pyplot as plt

files = {
    "DGAT": "./data/ARS_results/DGAT_Spearman.csv",
    "scLinear": "./data/ARS_results/scLinear_Spearman.csv",
    "CTP-net": "./data/ARS_results/CTP-net_Spearman.csv",
    "mRNA": "./data/ARS_results/mRNA_Spearman.csv",
    "sciPENN": "./data/ARS_results/sciPENN_Spearman.csv",
    "Seurat": "./data/ARS_results/Seurat_Spearman.csv"
}

dfs = []
for method, filepath in files.items():
    df = pd.read_csv(filepath)
    df["method"] = method
    dfs.append(df)

data = pd.concat(dfs, ignore_index=True)

# Rename columns for clarity
metric_cols = ["Protein_Spearman", "Protein_RMSE", "Cell_Spearman", "Cell_RMSE"]
for col in metric_cols:
    data[col] = pd.to_numeric(data[col], errors="coerce")

data.rename(columns={"Cell_Spearman": "Spot_Spearman", "Cell_RMSE": "Spot_RMSE"}, inplace=True)

def assign_scores(x, ascending):
    n = x.shape[0]
    ranks = x.rank(method="min", ascending=ascending)
    scores = n + 1 - ranks
    return scores

data["Spearman_protein_score"] = data.groupby("Sample")["Protein_Spearman"].transform(
    lambda x: assign_scores(x, ascending=False)
)
data["Spot_Spearman_score"] = data.groupby("Sample")["Spot_Spearman"].transform(
    lambda x: assign_scores(x, ascending=False)
)
data["Protein_RMSE_score"] = data.groupby("Sample")["Protein_RMSE"].transform(
    lambda x: assign_scores(x, ascending=True)
)
data["Spot_RMSE_score"] = data.groupby("Sample")["Spot_RMSE"].transform(
    lambda x: assign_scores(x, ascending=True)
)


agg_scores = data.groupby("method").agg({
    "Spearman_protein_score": "mean",
    "Spot_Spearman_score": "mean",
    "Protein_RMSE_score": "mean",
    "Spot_RMSE_score": "mean"
}).reset_index()

agg_scores["Protein_avg"] = (agg_scores["Spearman_protein_score"] + agg_scores["Protein_RMSE_score"]) / 2
agg_scores["Spot_avg"] = (agg_scores["Spot_Spearman_score"] + agg_scores["Spot_RMSE_score"]) / 2
agg_scores["Overall_ARScore"] = (agg_scores["Protein_avg"] + agg_scores["Spot_avg"]) / 2

print("Detailed Ranking Results (based on per-Sample scores):")
print(agg_scores[["method", "Spearman_protein_score", "Protein_RMSE_score", "Protein_avg",
                  "Spot_Spearman_score", "Spot_RMSE_score", "Spot_avg", "Overall_ARScore"]].to_string(index=False))


puor_colors = ["#8073ac", "#b35806","#b2abd2",
     "#e08214", "#fdb863",
    "#fee0b6"
]


n_methods = len(files)
colors = puor_colors[:n_methods]

def plot_barh(ax, methods, values, title, xlabel):
    sorted_idx = values.argsort()[::-1]
    values_sorted = values[sorted_idx]
    methods_sorted = methods[sorted_idx]
    ax.barh(methods_sorted, values_sorted, color=colors)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(xlabel)
    ax.invert_yaxis()  # Put the highest value on top

fig, axs = plt.subplots(1, 1, figsize=(4, 4))
fig.subplots_adjust(hspace=0.6)

# Plot Protein-level ranking scores (combining Spearman and RMSE for Protein)
plot_barh(
    #axs[0],
    axs,
    agg_scores["method"].values,
    agg_scores["Protein_avg"].values,
    "Protein-Level Ranking Score",
    "Ranking Score (Higher is Better)"
)

plt.tight_layout()
plt.savefig("Protein_Ranking_Score.pdf", dpi=600)
#plt.show()

fig, axs = plt.subplots(1, 1, figsize=(4, 4))
fig.subplots_adjust(hspace=0.6)

plot_barh(
    axs,
    agg_scores["method"].values,
    agg_scores["Spot_avg"].values,
    "Spot-Level", #Ranking Score (Spearman Corr & RMSE)",
    "Ranking Score (Higher is Better)"
)

plt.tight_layout()
plt.savefig("Figure_2D.pdf", dpi=600)
#plt.show()
