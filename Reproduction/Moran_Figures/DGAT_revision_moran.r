rm(list=ls())
options(stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("/Users/osmanbeyogluhu2/Documents/Manuscripts/DGAT/manuscript/NatComm/revision/update/")
moran_df <- read.csv("Moran_DGAT_mRNA.csv", row.names = 1)

dat_tnbc <- moran_df %>% filter(tissue == "TNBC")

#------------------------------------------------------------
# Mapping markers → original modules
#------------------------------------------------------------
dat_tnbc <- dat_tnbc %>%
  mutate(module = case_when(
    marker %in% c("CD14","CD68","CD163","FCGR3A","ITGAM","ITGAX","CEACAM8") ~ "Myeloid / APC",
    marker %in% c("CD3E","CD4","CD8A","CCR7","PDCD1","CD27")                 ~ "T cell",
    marker %in% c("CD19","PAX5","MS4A1","CD40","CR2","CXCR5")               ~ "B cell / TLS",
    marker %in% c("EPCAM","KRT5","SDC1")                                    ~ "Basal epithelial",
    marker %in% c("ACTA2","VIM","PECAM1")                                   ~ "Stromal / endothelial",
    marker %in% c("CD274")                                                  ~ "Immune checkpoint",
    marker %in% c("PCNA")                                                   ~ "Proliferation",
    TRUE                                                                    ~ "Other"
  ))

#------------------------------------------------------------
# Merge categories for simpler legend
#------------------------------------------------------------
dat_tnbc <- dat_tnbc %>%
  mutate(module_merged = case_when(
    module %in% c("T cell", "Immune checkpoint")            ~ "T cell / checkpoint",
    module %in% c("Proliferation", "Other")                 ~ "Other",
    TRUE                                                    ~ module
  ))

#------------------------------------------------------------
# Color palette for merged groups
#------------------------------------------------------------
tnbc_palette_merged <- c(
  "Myeloid / APC"          = "#E66101",
  "T cell / checkpoint"    = "#1F78B4",
  "B cell / TLS"           = "#33A02C",
  "Basal epithelial"       = "#6A3D9A",
  "Stromal / endothelial"  = "#B15928",
  "Other"                  = "grey60"
)

#------------------------------------------------------------
# FINAL PLOT with bottom legend (2 rows)
#------------------------------------------------------------
p_tnbc <- ggplot(dat_tnbc, aes(x = moran_mrna, y = moran_protein)) +
  
  geom_abline(slope = 1, intercept = 0,
              aes(linetype = "Identity (y = x)"),
              color = "grey50", linewidth = 0.5) +
  
  geom_point(aes(color = module_merged), size = 3, alpha = 0.9) +
  
  geom_text_repel(
    aes(label = marker, color = module_merged),
    show.legend = FALSE,
    size = 3.5,
    max.overlaps = Inf,
    point.padding = 0.2,
    box.padding = 0.3,
    segment.size = 0.2,
    seed = 123
  ) +
  
  scale_x_continuous("mRNA Moran's I",
                     limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  
  scale_y_continuous("DGAT Moran's I",
                     limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     expand = c(0, 0)) +
  
  coord_equal() +
  
  scale_color_manual(values = tnbc_palette_merged,
                     name = "Marker type") +
  
  scale_linetype_manual(values = c("Identity (y = x)" = "dashed"),
                        name = "") +
  
  #ggtitle("Spatial autocorrelation (Moran’s I)") +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 15),
    #axis.title   = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 12),     # <-- legend item text size
    legend.spacing.x = unit(3, "mm")
  ) +
  
  guides(
    color = guide_legend(title = "", nrow = 3, byrow = TRUE),
    linetype = guide_legend(order = 2)
  ) 
p_tnbc


#-----------------------------
# 6. Save publication-quality files
#-----------------------------
ggsave("TNBC_MoransI_DGAT_vs_mRNA.pdf",
       p_tnbc, width = 6, height = 6, units = "in")

ggsave("TNBC_MoransI_DGAT_vs_mRNA.png",
       p_tnbc, width = 6, height = 6, units = "in", dpi = 600)
