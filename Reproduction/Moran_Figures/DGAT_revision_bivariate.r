rm(list=ls())
options(stringsAsFactors = FALSE)
#############################################
## DGAT TNBC: Bivariate Moran + Network
#############################################
library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(igraph)
library(ggraph)
library(patchwork)

## -----------------------------
## 0. Setup & data loading
## -----------------------------

# Set your working directory (adjust as needed)
setwd("/Users/osmanbeyogluhu2/Documents/Manuscripts/DGAT/manuscript/NatComm/revision/update/Bivariate/")

# Read DGAT bivariate Moran's I results for TNBC
#df_protein <- read_csv("DGAT_Bivariate_Moran_Unique/Bivariate_Unique_TNBC.csv")
#df_protein <- read_csv("DGAT_Bivariate_Moran_Unique/Bivariate_Unique_iBreast.csv")
df_protein <- read_csv("DGAT_Bivariate_Moran_Unique/Bivariate_Unique_Melanoma.csv")
#df_protein <- read_csv("DGAT_Bivariate_Moran_Unique/Bivariate_Unique_Prostate.csv")

# Optional: Define marker types for coloring the network
marker_info <- tibble::tibble(
  marker = c("ACTA2","BCL2","CCR7","CD14","CD163","CD19","CD27","CD274","CD3E",
             "CD4","CD40","CD68","CD8A","CEACAM8","CR2","CXCR5","EPCAM","FCGR3A",
             "ITGAM","ITGAX","KRT5","MS4A1","PAX5","PCNA","PDCD1","PECAM1","SDC1","VIM"),
  type   = c("Stromal / endothelial",  # ACTA2
             "Basal epithelial",       # BCL2
             "T cell",                 # CCR7
             "Myeloid / APC",          # CD14
             "Myeloid / APC",          # CD163
             "B cell / TLS",           # CD19
             "T cell",                 # CD27
             "Immune checkpoint",      # CD274
             "T cell",                 # CD3E
             "T cell",                 # CD4
             "B cell / TLS",           # CD40
             "Myeloid / APC",          # CD68
             "T cell",                 # CD8A
             "Myeloid / APC",          # CEACAM8
             "B cell / TLS",           # CR2
             "B cell / TLS",           # CXCR5
             "Basal epithelial",       # EPCAM
             "Myeloid / APC",          # FCGR3A
             "Myeloid / APC",          # ITGAM
             "Myeloid / APC",          # ITGAX
             "Basal epithelial",       # KRT5
             "B cell / TLS",           # MS4A1
             "B cell / TLS",           # PAX5
             "Proliferation",          # PCNA
             "T cell",                 # PDCD1
             "Stromal / endothelial",  # PECAM1
             "Basal epithelial",       # SDC1
             "Stromal / endothelial"   # VIM
  )
)

type_cols <- c(
  "B cell / TLS"          = "#33A02C",
  "Basal epithelial"      = "#6A3D9A",
  "Immune checkpoint"     = "#D73027",
  "Myeloid / APC"         = "#E66101",
  "Proliferation"         = "#FFD92F",
  "Stromal / endothelial" = "#B15928",
  "T cell"                = "#1F78B4",
  "Other"                 = "grey70"
)

## -----------------------------
## 1. Build symmetric matrix
## -----------------------------

protein_list <- sort(unique(c(df_protein$Marker_A, df_protein$Marker_B)))

mat <- matrix(NA_real_,
              nrow = length(protein_list),
              ncol = length(protein_list),
              dimnames = list(protein_list, protein_list))

for (i in seq_len(nrow(df_protein))) {
  a <- df_protein$Marker_A[i]
  b <- df_protein$Marker_B[i]
  val <- df_protein$Bivariate_Moran_I[i]
  mat[a, b] <- val
  mat[b, a] <- val
}

# Cluster markers for ordering
mat_for_dist <- mat
mat_for_dist[is.na(mat_for_dist)] <- 0
hc <- hclust(dist(mat_for_dist), method = "average")
ordered <- rownames(mat_for_dist)[hc$order]
mat_ordered <- mat[ordered, ordered]

## -----------------------------
## 2. Lower-triangle heatmap
## -----------------------------

## -----------------------------
## Heatmap rotated 270 degrees
## -----------------------------

df_long <- melt(mat_ordered, na.rm = FALSE)
colnames(df_long) <- c("Marker_B", "Marker_A", "Moran_I")

df_tri <- df_long %>%
  mutate(
    Marker_A = factor(Marker_A, levels = ordered),
    Marker_B = factor(Marker_B, levels = ordered),
    i = as.numeric(Marker_A),
    j = as.numeric(Marker_B)
  ) %>%
  filter(j >= i) %>%       # lower triangle before rotation
  select(-i, -j)

threshold <- 0.70
df_high <- df_tri %>% filter(!is.na(Moran_I) & Moran_I >= threshold)

heatmap_plot <- ggplot(df_tri,
                       aes(x = Marker_B,     # swap axes
                           y = Marker_A,     # swap axes
                           fill = Moran_I)) +

  geom_tile(color = "grey90", linewidth = 0.2) +

  geom_text(
    data = df_high,
    aes(label = sprintf("%.2f", Moran_I)),
    color = "black",
    size = 2.5,
    fontface = "bold",
    angle = 0,                 # rotate labels as in your PNG
    vjust = 0.5,
    hjust = 0.5
  ) +

  scale_fill_gradientn(
    colours  = c("#4575b4", "white", "#f46d43"),
    values   = c(0, 0.5, 1),
    limits   = c(0, 1),
    name     = "Bivariate Moran’s",
    na.value = "white"
  ) +

  coord_equal(expand = FALSE) +

  scale_x_discrete(limits = rev(ordered)) +   # *** rotate 270° ***
  scale_y_discrete(limits = ordered) +

  #labs(title = "Bivariate Moran’s") +

  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),

    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 8),

    axis.line.x  = element_line(color = "black"),
    axis.line.y  = element_line(color = "black"),

    legend.position = "right",
    legend.title    = element_text(face = "bold")
  )

heatmap_plot

library(grid)

# build the legend as a separate grob
legend_grob <- cowplot::get_legend(
  heatmap_plot +
    guides(fill = guide_colorbar(
      title = "Bivariate Moran’s I",
      barheight = unit(40, "mm"),
      barwidth  = unit(6, "mm")
    )) +
    theme(legend.position = "right",
          legend.title = element_text(face = "bold"))
)

# remove default legend from heatmap
heatmap_noleg <- heatmap_plot + theme(legend.position = "none")

# place legend manually into the empty space
heatmap_final <- heatmap_noleg +
  annotation_custom(
    grob = legend_grob,
    xmin = length(ordered)*0.55,   # adjust horizontal position
    xmax = length(ordered)*0.95,   # adjust width
    ymin = length(ordered)*0.55,   # adjust vertical position
    ymax = length(ordered)*0.95
  )
heatmap_final
## -----------------------------
## 3. Network from strong pairs
## -----------------------------
# --- build graph as before ---
df_protein_unique <- df_protein %>%
  rowwise() %>%
  mutate(pair = paste(sort(c(Marker_A, Marker_B)), collapse = "_")) %>%
  ungroup() %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

df_final <- df_protein_unique %>%
  filter(Bivariate_Moran_I >= threshold)

g <- graph_from_data_frame(
  df_final %>% select(Marker_A, Marker_B, Bivariate_Moran_I),
  directed = FALSE
)

E(g)$weight <- df_final$Bivariate_Moran_I

marker_type_vec <- marker_info$type
names(marker_type_vec) <- marker_info$marker
V(g)$type <- marker_type_vec[V(g)$name]
V(g)$type[is.na(V(g)$type)] <- "Other"

set.seed(123)
layout_fr <- create_layout(g, layout = "fr")

## make layout a bit narrower
layout_fr$x <- layout_fr$x * 0.5

network <- ggraph(layout_fr) +
  geom_edge_link(aes(width = weight, colour = weight),
                 alpha = 0.8, show.legend = TRUE) +
  geom_node_point(aes(colour = type), size = 4, show.legend = TRUE) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_width(range = c(0.5, 2.5), guide = "none") +
  # edge colourbar
  scale_edge_colour_gradientn(
    colours = c("#4575b4", "#f46d43"),
    limits  = c(threshold, 1),
    name    = "Bivariate Moran’s I",
    guide   = guide_colorbar(order = 1)
  ) +
  # node colours
  scale_colour_manual(
    values = type_cols,
   # name   = "Marker type",
    guide  = guide_legend(order = 2)
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.margin  = margin(5, 5, 5, 5),
    legend.position = "bottom",
    legend.box      = "vertical",      # colourbar above, types below
    legend.title    = element_text(face = "bold", size = 10),
    legend.text     = element_text(size = 9),
    legend.key.size = unit(3, "mm")
  )

network

# example narrow save
# ggsave("network_narrow.pdf", network, width = 2.8, height = 5.5, units = "in", dpi = 300)

## -----------------------------
## 4. Combine heatmap + network
## -----------------------------

heatmap_clean <- heatmap_final +
  theme(plot.margin = margin(5, 5, 5, 5))

network_clean <- network +
  theme(plot.margin = margin(5, 5, 5, 5))

combined_tnbc_assoc <- heatmap_clean + network_clean +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(
    tag_levels = list(c("f", "g")),
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12),
      plot.margin = margin(5, 5, 5, 5),
      text = element_text(size = 11)
    )
  )

# Show combined figure
combined_tnbc_assoc

# Save high-res version for manuscript
ggsave(
  "TNBC_bivariate_Moran_and_network.png",
  combined_tnbc_assoc,
  width = 8,
  height = 4,
  dpi = 600
)



# rm(list=ls())
# options(stringsAsFactors = FALSE)
# # Load necessary libraries
# library(readr)
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(patchwork)
# library(stringr)
# 
# # --- 1. Data Loading and Initial Setup ---
# setwd("/Users/osmanbeyogluhu2/Documents/Manuscripts/DGAT/manuscript/NatComm/revision/update/Bivariate/")
# # Load the DGAT Predicted Protein data (use uploaded file name)
# df_protein <- read_csv("DGAT_Bivariate_Moran_Unique/Bivariate_Unique_TNBC.csv") %>%
#   mutate(Dataset = "Predicted Protein (DGAT)")
# 
# # Load the mRNA/Comparison data (use uploaded file name)
# df_mrna <- read_csv("mRNA_Bivariate_Moran_Unique/Bivariate_Unique_TNBC.csv") %>%
#   mutate(Dataset = "mRNA Baseline")
# 
# # Combine for systemic analysis (Panel A)
# df_combined <- bind_rows(df_protein, df_mrna)
# 
# 
# # library(ggplot2)
# # library(dplyr)
# # library(reshape2)
# # 
# # # --- 1. Set threshold for strong correlations ---
# # threshold <- 0.7
# # 
# # # --- 2. Prepare square matrix of proteins ---
# # protein_list <- unique(c(df_protein$Marker_A, df_protein$Marker_B))
# # mat <- matrix(NA, nrow=length(protein_list), ncol=length(protein_list),
# #               dimnames=list(protein_list, protein_list))
# # 
# # # Fill matrix with Moran's I values
# # for(i in 1:nrow(df_protein)){
# #   a <- df_protein$Marker_A[i]
# #   b <- df_protein$Marker_B[i]
# #   mat[a,b] <- df_protein$Bivariate_Moran_I[i]
# #   mat[b,a] <- df_protein$Bivariate_Moran_I[i]  # symmetric
# # }
# # 
# # --- 3. Hierarchical clustering ---
# row_order <- hclust(dist(mat))$order
# col_order <- hclust(dist(t(mat)))$order
# mat_ordered <- mat[row_order, col_order]
# 
# # --- 4. Convert to long format ---
# df_long <- melt(mat_ordered, na.rm=FALSE)
# colnames(df_long) <- c("Marker_A", "Marker_B", "Bivariate_Moran_I")
# 
# # --- 5. Keep only upper triangle ---
# df_long <- df_long %>%
#   filter(as.numeric(factor(Marker_A)) <= as.numeric(factor(Marker_B)))
# 
# # --- 6. Identify strong pairs and top 3 ---
# strong_pairs <- df_long %>% filter(Bivariate_Moran_I >= threshold)
# top3_pairs <- strong_pairs %>% arrange(desc(Bivariate_Moran_I)) %>% slice(1:3)
# 
# # --- 7. Plot heatmap ---
# heatmap_DGAT <- ggplot(df_long, aes(x=Marker_A, y=Marker_B, fill=Bivariate_Moran_I)) +
#   geom_tile(color="white") +
#   geom_tile(data=strong_pairs, aes(x=Marker_A, y=Marker_B),
#             fill="black", alpha=0.5, width=1, height=1) +
#   geom_text(data=strong_pairs, aes(label=round(Bivariate_Moran_I,2)),
#             color="white", size=3) +
#   geom_text(data=top3_pairs, aes(label=round(Bivariate_Moran_I,2)),
#             color="white", fontface="bold", size=3.5) +
#   scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0, na.value="white") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=45, hjust=1),
#         axis.text.y = element_text(angle=0)) +
#   labs(title=paste0("DGAT Predicted Protein: Moran's I ≥ ", threshold),
#        fill="Moran I")
# 
# heatmap_DGAT
# #################
# library(dplyr)
# library(igraph)
# library(ggraph)
# 
# # Remove duplicate protein pairs
# df_protein_unique <- df_protein %>%
#   rowwise() %>%
#   mutate(pair = paste(sort(c(Marker_A, Marker_B)), collapse = "_")) %>%
#   ungroup() %>%
#   distinct(pair, .keep_all = TRUE) %>%
#   select(-pair)
# 
# # Build final network for cutoff 0.7
# df_final <- df_protein_unique %>%
#   filter(Bivariate_Moran_I >= 0.7)
# 
# # Create igraph object
# g <- graph_from_data_frame(
#   df_final %>% select(Marker_A, Marker_B), 
#   directed = FALSE
# )
# 
# # Assign edge weights BEFORE simplifying
# E(g)$weight <- df_final$Bivariate_Moran_I
# 
# # Simplify network but KEEP edge attributes
# g <- simplify(g, edge.attr.comb = list(weight = "mean"))
# 
# # Community detection
# comm <- cluster_fast_greedy(g, weights = E(g)$weight)
# V(g)$community <- membership(comm)
# 
# # Plot
# network <- ggraph(g, layout = "fr") +
#   geom_edge_link(aes(width = weight, color = weight), alpha = 0.85) +
#   geom_node_point(aes(color = factor(community)), size = 5) +
#   geom_node_text(aes(label = name), repel = TRUE, size = 4) +
#   scale_edge_color_distiller(palette = "RdBu", limits = c(-1, 1), direction = 1) +
#   theme_void() +
#   ggtitle("Spatial Protein Association Network (DGAT, TNBC)")
# 
# heatmap_DGAT + network
