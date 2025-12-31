# ---- Load packages ----
library(tidyverse)
library(ggplot2)

# ---- Step 1: Load data ----
setwd("/Users/osmanbeyogluhu2/Documents/Manuscripts/DGAT/manuscript/NatComm/revision/")
df <- read.csv("DGAT_mRNA_samplewise.csv",
               header = TRUE, check.names = FALSE)
rownames(df) <- df[, 1]
df <- df[, -1, drop = FALSE]

# ---- Step 2: Define marker-to-group mapping ----
marker_groups <- list(
  T_cell        = c("CD3E","CD4","CD8A","PDCD1","CD27","CCR7","CXCR5"),
  B_cell        = c("CD19","MS4A1","PAX5","CR2","BCL2"),
  Myeloid       = c("CD14","CD68","ITGAM","ITGAX","FCGR3A","CEACAM8","CD163"),
  Antigen_presenting = c("HLA_DRA","CD40"),
  Epithelial    = c("EPCAM","KRT5"),
  Endothelial   = c("PECAM1"),
  Stromal       = c("ACTA2","VIM","SDC1"),
  Proliferative = c("PCNA"),
  Checkpoint    = c("CD274"),
  Leukocyte     = c("PTPRC_1","PTPRC_2")
)

# Build mapping table
marker_to_group <- enframe(marker_groups, name = "group", value = "marker") %>%
  unnest(marker) %>%
  filter(marker %in% colnames(df))

# Optional: collapse replicate CD45 features (PTPRC_1/PTPRC_2) to a single "PTPRC"
# If you prefer to keep them separate, set COLLAPSE_PTPRC <- FALSE
COLLAPSE_PTPRC <- FALSE
if (COLLAPSE_PTPRC) {
  # If both exist, replace them with their average and rename to "PTPRC"
  ptprc_feats <- intersect(c("PTPRC_1","PTPRC_2"), colnames(df))
  if (length(ptprc_feats) >= 1) {
    df$PTPRC <- rowMeans(df[, ptprc_feats, drop = FALSE], na.rm = TRUE)
    df <- df %>% select(-all_of(ptprc_feats))
    marker_to_group <- marker_to_group %>%
      filter(!marker %in% ptprc_feats) %>%
      bind_rows(tibble(group = "Leukocyte", marker = "PTPRC"))
  }
}

# ---- Step 3: Reshape data ----
# NOTE: This assumes rownames are something like "Tonsil1_DGAT" or "Tonsil1_mRNA"
# We parse Model from the Sample string and remove it to get Tissue.
df_long <- df %>%
  rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "marker", values_to = "correlation") %>%
  mutate(
    Model = case_when(
      grepl("DGAT", Sample, ignore.case = TRUE) ~ "DGAT",
      grepl("mRNA", Sample, ignore.case = TRUE) ~ "mRNA",
      TRUE ~ NA_character_
    ),
    Tissue = Sample %>%
      gsub("(?i)_?DGAT.*$", "", ., perl = TRUE) %>%  # remove DGAT suffix if present
      gsub("(?i)_?mRNA.*$", "", ., perl = TRUE) %>%  # remove mRNA suffix if present
      gsub("_.*$", "", .),                           # keep first token as tissue if needed
    group = marker_to_group$group[match(marker, marker_to_group$marker)]
  ) %>%
  filter(!is.na(group), !is.na(Model))

# Set a consistent group order (matches your manuscript language)
group_order <- c("T_cell","B_cell","Myeloid","Antigen_presenting",
                 "Epithelial","Endothelial","Stromal","Proliferative",
                 "Checkpoint","Leukocyte")
df_long <- df_long %>%
  mutate(group = factor(group, levels = group_order))

# ---- Step 4: Boxplot ----
ggplot(df_long, aes(x = group, y = correlation, fill = Model)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85,
               position = position_dodge(width = 0.8), width = 0.7) +
  geom_jitter(aes(color = Model),
              alpha = 0.55, size = 1.2,
              position = position_jitterdodge(jitter.width = 0.18, dodge.width = 0.8)) +
  facet_wrap(~ Tissue, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("DGAT" = "#7B68EE", "mRNA" = "#FFA500")) +
  scale_color_manual(values = c("DGAT" = "#7B68EE", "mRNA" = "#FFA500")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "DGAT-predicted protein vs mRNA cell type marker groups across tissue types",
    x = " ",
    y = "Spearman Correlation",
    fill = "Model",
    color = "Model"
  )



