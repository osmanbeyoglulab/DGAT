# ============================================================
# Supplementary Table X (ALL proteins)
# - Starts from your input: sample_protein_top20_gene_correlations-2.csv
# - Keeps only gene–protein pairs with mean |r| >= r_cut across >= min_tissues
# - Reports up to top_n genes per protein AFTER filtering
# ============================================================

library(tidyverse)

# -----------------------------
# INPUT
# -----------------------------
setwd("/Users/osmanbeyogluhu2/Documents/Manuscripts/DGAT/manuscript/NatComm/revision/")
df <- read.csv("sample_protein_top20_gene_correlations-2.csv", check.names = FALSE)

# -----------------------------
# REQUIRED: training/eval mapping (use your existing set_map)
# Example:
# set_map <- c("Tonsil"="Training", ... , "CID44971"="Evaluation", ...)
# names(set_map) must match df$Sample values
# -----------------------------

df_tbl <- df %>%
  mutate(
    Correlation = as.numeric(Correlation),
    Dataset = set_map[Sample],
    Dataset = factor(Dataset, levels = c("Training", "Evaluation"))
  ) %>%
  filter(!is.na(Dataset), !is.na(Correlation))

# -----------------------------
# PARAMETERS (EDIT IF NEEDED)
# -----------------------------
r_cut       <- 0.2   # correlation strength cutoff
min_tissues <- 4      # require consistency across tissues
top_n       <- 20     # keep at most 20 genes per protein after filtering

# -----------------------------
# 1) Compute pair-level stats across tissues
# -----------------------------
pair_stats <- df_tbl %>%
  group_by(Protein, mRNA) %>%
  summarise(
    n_tissues_total = n_distinct(Sample),
    mean_abs_r      = mean(abs(Correlation), na.rm = TRUE),
    mean_r          = mean(Correlation, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 2) Filter by mean |r| and tissue consistency
# -----------------------------
filtered_pairs <- pair_stats %>%
  filter(n_tissues_total >= min_tissues, mean_abs_r >= r_cut)

# -----------------------------
# 3) Select up to top_n genes per protein (after filtering)
#    Ranked by: consistency (n_tissues_total), then strength (mean_abs_r)
# -----------------------------
top_pairs <- filtered_pairs %>%
  group_by(Protein) %>%
  arrange(desc(n_tissues_total), desc(mean_abs_r)) %>%
  slice_head(n = top_n) %>%
  ungroup()

# -----------------------------
# 4) Summaries by dataset (training vs evaluation)
# -----------------------------
median_by_dataset <- df_tbl %>%
  semi_join(top_pairs, by = c("Protein", "mRNA")) %>%
  group_by(Protein, Dataset) %>%
  summarise(
    Median_r = median(Correlation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Dataset,
    values_from = Median_r,
    names_prefix = "Median_r_"
  )

tissue_by_dataset <- df_tbl %>%
  semi_join(top_pairs, by = c("Protein", "mRNA")) %>%
  group_by(Protein, Dataset) %>%
  summarise(
    n_tissues = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Dataset,
    values_from = n_tissues,
    names_prefix = "n_tissues_"
  )

# -----------------------------
# 5) Collapse gene lists per protein
# -----------------------------
gene_list <- top_pairs %>%
  group_by(Protein) %>%
  summarise(
    n_genes = n(),
    Top_correlated_genes = paste(mRNA, collapse = ", "),
    .groups = "drop"
  )

# -----------------------------
# 6) Final table
# -----------------------------
supp_table <- gene_list %>%
  left_join(median_by_dataset, by = "Protein") %>%
  left_join(tissue_by_dataset, by = "Protein") %>%
  arrange(Protein)

# -----------------------------
# 7) Save
# -----------------------------
out_file <- sprintf(
  "Supplementary_Table_X_Filtered_top%d_meanAbsR_ge_%.2f_minT%d.csv",
  top_n, r_cut, min_tissues
)

write.csv(supp_table, out_file, row.names = FALSE)
message("Saved: ", out_file)
