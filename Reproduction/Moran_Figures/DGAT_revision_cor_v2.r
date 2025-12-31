rm(list=ls())
options(stringsAsFactors = FALSE)

library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# -----------------------------
# helper: normalize dash characters so tissue/sample names match
# -----------------------------
norm_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", x, perl = TRUE)
  x
}

# -----------------------------
# 0. Load data
# -----------------------------
setwd("/Users/osmanbeyogluhu2/Documents/Manuscripts/DGAT/manuscript/NatComm/revision/")
df <- read.csv("sample_protein_top20_gene_correlations-2.csv", check.names = FALSE)

setwd("/Users/osmanbeyogluhu2/Documents/Manuscripts/DGAT/manuscript/NatComm/revision/update/")
moran_df <- read.csv("Moran_DGAT_mRNA.csv", check.names = FALSE)

# -----------------------------
# 1. Parameters (set once)
# -----------------------------
corr_thr    <- 0.2
min_tissues <- 4
top_genes   <- 20
lim_fixed   <- 0.7           # fixed max for 0->lim color scale (NULL for per-protein autoscale)
na_col      <- "grey95"

# -----------------------------
# 2. Sample order, display labels, and train/eval mapping
# -----------------------------
sample_order <- norm_name(c(
  "Tonsil",
  "Tonsil_AddOns",
  "PBC-PR_6835-5A",
  "PBC_PR_6837",
  "Glioblastoma",
  "Breast",
  "CID44971",
  "LN",
  "Melanoma",
  "Prostate",
  "iBreast_test"
))

# display labels (for figures only)
sample_label_map <- c(
  "Tonsil"          = "Tonsil1",
  "Tonsil_AddOns"   = "Tonsil2",
  "PBC-PR_6835-5A"  = "Meso1",
  "PBC_PR_6837"     = "Meso2",
  "Glioblastoma"    = "GBM",
  "Breast"          = "BC",
  "CID44971"        = "TNBC",
  "LN"              = "LN",
  "Melanoma"        = "Melanoma",
  "Prostate"        = "PC",
  "iBreast_test"    = "BC_test"
)
names(sample_label_map) <- norm_name(names(sample_label_map))

set_map <- c(
  "Tonsil"         = "Training",
  "Tonsil_AddOns"  = "Training",
  "PBC-PR_6835-5A" = "Training",
  "PBC_PR_6837"    = "Training",
  "Glioblastoma"   = "Training",
  "Breast"         = "Training",
  "CID44971"       = "Evaluation",
  "LN"             = "Evaluation",
  "Melanoma"       = "Evaluation",
  "Prostate"       = "Evaluation",
  "iBreast_test"   = "Evaluation"
)
names(set_map) <- norm_name(names(set_map))

df <- df %>%
  mutate(
    Sample = norm_name(Sample),
    Sample = factor(Sample, levels = sample_order),
    Protein = as.character(Protein),
    mRNA = as.character(mRNA),
    Correlation = as.numeric(Correlation)
  )

# -----------------------------
# 3. Moran lookup (Protein x Sample -> moran_protein)
# -----------------------------
cn <- colnames(moran_df)
bad <- is.na(cn) | cn == ""
if (any(bad)) {
  cn[bad] <- paste0("Xbad", seq_len(sum(bad)))
  colnames(moran_df) <- cn
}
colnames(moran_df) <- trimws(colnames(moran_df))

# expected columns: tissue, marker, moran_protein, moran_mrna
moran_df <- moran_df %>%
  rename(
    tissue        = tissue,
    marker        = marker,
    moran_protein = moran_protein,
    moran_mrna    = moran_mrna
  )

tissue_map <- c("Lymph Node" = "LN")
names(tissue_map) <- norm_name(names(tissue_map))
tissue_map <- setNames(norm_name(tissue_map), names(tissue_map))

moran_lookup <- moran_df %>%
  mutate(
    Sample  = norm_name(tissue),
    Sample  = dplyr::recode(Sample, !!!tissue_map, .default = Sample),
    Sample  = factor(Sample, levels = sample_order),
    Protein = as.character(marker),
    moran_protein = as.numeric(moran_protein)
  ) %>%
  filter(!is.na(Sample)) %>%
  select(Protein, Sample, moran_protein) %>%
  distinct()

# -----------------------------
# 4. Select consistent genes per protein
# -----------------------------
consistent_genes <- df %>%
  filter(abs(Correlation) >= corr_thr) %>%
  group_by(Protein, mRNA) %>%
  summarise(
    n_tissue  = n_distinct(Sample),
    mean_corr = mean(Correlation, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  group_by(Protein) %>%
  filter(n_tissue >= min_tissues) %>%
  arrange(desc(n_tissue), desc(abs(mean_corr))) %>%
  slice_head(n = top_genes) %>%
  ungroup()

df_consistent <- df %>%
  semi_join(consistent_genes, by = c("Protein", "mRNA"))

# -----------------------------
# 5. Heatmap function (positive-only; NA-aware; split Train/Eval; updated labels)
# -----------------------------
make_heatmap <- function(protein_name, show_stats = TRUE, split_train_eval = TRUE) {
  
  df_prot <- df_consistent %>% filter(Protein == protein_name)
  if (nrow(df_prot) == 0) return(NULL)
  
  df_prot_agg <- df_prot %>%
    group_by(mRNA, Sample) %>%
    summarise(Correlation = mean(as.numeric(Correlation), na.rm = TRUE), .groups = "drop") %>%
    mutate(Sample = as.character(Sample))
  
  mat_df <- df_prot_agg %>%
    select(mRNA, Sample, Correlation) %>%
    pivot_wider(names_from = Sample, values_from = Correlation, values_fill = NA_real_) %>%
    as.data.frame(check.names = FALSE)
  
  rownames(mat_df) <- mat_df$mRNA
  mat_df$mRNA <- NULL
  mat <- as.matrix(mat_df)
  
  # ensure all sample_order columns exist and ordered
  missing_cols <- setdiff(sample_order, colnames(mat))
  if (length(missing_cols) > 0) {
    add_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = length(missing_cols),
                      dimnames = list(rownames(mat), missing_cols))
    mat <- cbind(mat, add_mat)
  }
  mat <- mat[, sample_order, drop = FALSE]
  
  # dataset vec (internal names, ordered)
  dataset_vec_internal <- unname(set_map[sample_order])
  names(dataset_vec_internal) <- sample_order
  
  # order rows by mean correlation in training
  train_cols <- sample_order[dataset_vec_internal == "Training"]
  row_score <- rowMeans(mat[, train_cols, drop = FALSE], na.rm = TRUE)
  ord <- order(row_score, decreasing = TRUE, na.last = TRUE)
  mat <- mat[ord, , drop = FALSE]
  
  # convert column names to display labels
  colnames(mat) <- sample_label_map[colnames(mat)]
  display_order <- sample_label_map[sample_order]
  mat <- mat[, display_order, drop = FALSE]
  
  # dataset vec in display names (for split/annotation)
  dataset_vec <- unname(set_map[sample_order])
  names(dataset_vec) <- display_order
  
  # optional stats line
  stats_label <- NULL
  if (show_stats) {
    train_disp <- display_order[dataset_vec == "Training"]
    eval_disp  <- display_order[dataset_vec == "Evaluation"]
    
    n_genes <- nrow(mat)
    n_eval_present <- sum(rowSums(!is.na(mat[, eval_disp, drop = FALSE])) > 0)
    
    med_train <- suppressWarnings(median(as.vector(mat[, train_disp, drop = FALSE]), na.rm = TRUE))
    med_eval  <- suppressWarnings(median(as.vector(mat[, eval_disp,  drop = FALSE]), na.rm = TRUE))
    
    stats_label <- sprintf(
      "median r: train=%.2f, eval=%.2f",
      #n_genes, n_eval_present, 100 * n_eval_present / max(n_genes, 1),
      med_train, med_eval
    )
  }
  
  # color scale (positive-only)
  lim <- lim_fixed
  if (is.null(lim)) {
    lim <- suppressWarnings(max(mat, na.rm = TRUE))
    if (!is.finite(lim) || lim <= 0) lim <- 1
  }
  col_corr <- colorRamp2(c(0, lim), c("#F7F7F7", "#B2182B"))
  
  # annotation: clarify paired vs predicted (legend labels)
  dataset_label <- ifelse(
    dataset_vec == "Training",
    "Training (paired mRNA+protein)",
    "Evaluation (mRNA-only; DGAT-pred protein)"
  )
  names(dataset_label) <- names(dataset_vec)
  
  ha <- HeatmapAnnotation(
    dataset = dataset_label,
    col = list(dataset = c(
      "Training (paired mRNA+protein)" = "#6BAED6",
      "Evaluation (mRNA-only; DGAT-pred protein)" = "#FD8D3C"
    )),
    #annotation_name_side = "left",
    annotation_name = NULL
  )
  
  Heatmap(
    mat,
    name = "Correlation (r)",
    col = col_corr,
    na_col = na_col,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = if (split_train_eval) factor(dataset_vec, levels = c("Training","Evaluation")) else NULL,
    top_annotation = ha,
    column_title = protein_name,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 9)
  )
}

# -----------------------------
# 6. Save one heatmap per protein
# -----------------------------
setwd("/Users/osmanbeyogluhu2/Documents/Manuscripts/DGAT/manuscript/NatComm/revision/update/")

proteins_to_plot <- sort(unique(df$Protein))

for (p in proteins_to_plot) {
  ht <- make_heatmap(p, show_stats = TRUE, split_train_eval = TRUE)
  if (!is.null(ht)) {
    message("Saving heatmap for ", p)
    pdf(paste0("Heatmap_", p, ".pdf"), width = 5.5, height = 4)
    draw(ht, merge_legend = TRUE)
    dev.off()
  } else {
    message("Skipping ", p, " (no genes after filtering).")
  }
}
########################################################################################################

library(tidyverse)

# -----------------------------
# Parameters
# -----------------------------
top_n_genes <- 20   # how many top genes per protein to report

# -----------------------------
# Add dataset label
# -----------------------------
df_tbl <- df %>%
  mutate(
    Dataset = set_map[Sample],
    Dataset = factor(Dataset, levels = c("Training", "Evaluation"))
  )

# -----------------------------
# Compute per-protein, per-gene summaries
# -----------------------------
gene_summary <- df_tbl %>%
  group_by(Protein, mRNA, Dataset) %>%
  summarise(
    median_r = median(Correlation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Dataset,
    values_from = median_r,
    names_prefix = "median_r_"
  )

# -----------------------------
# Rank genes by TRAINING correlation
# -----------------------------
top_genes <- gene_summary %>%
  arrange(Protein, desc(median_r_Training)) %>%
  group_by(Protein) %>%
  slice_head(n = top_n_genes) %>%
  ungroup()

# -----------------------------
# Count evaluation presence
# -----------------------------
eval_presence <- df_tbl %>%
  filter(Dataset == "Evaluation") %>%
  group_by(Protein, mRNA) %>%
  summarise(
    n_eval_samples = sum(!is.na(Correlation)),
    .groups = "drop"
  )

# -----------------------------
# Combine
# -----------------------------
top_genes_final <- top_genes %>%
  left_join(eval_presence, by = c("Protein", "mRNA")) %>%
  mutate(
    n_eval_samples = replace_na(n_eval_samples, 0)
  )

# -----------------------------
# Collapse genes into one row per protein
# -----------------------------
supp_table <- top_genes_final %>%
  group_by(Protein) %>%
  summarise(
    Top_correlated_genes = paste(mRNA, collapse = ", "),
    Median_r_Training   = round(median(median_r_Training, na.rm = TRUE), 3),
    Median_r_Evaluation = round(median(median_r_Evaluation, na.rm = TRUE), 3),
    .groups = "drop"
  )

# -----------------------------
# Save table
# -----------------------------
write.csv(
  supp_table,
  "Supplementary_Table_mRNA_Protein_Correlations.csv",
  row.names = FALSE
)



