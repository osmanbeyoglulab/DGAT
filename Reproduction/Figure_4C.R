library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(tibble)

FONTSIZE = 18
RELSIZE = 0.8
var = 'diff_zscore_global'

data = read.csv('./data/dotplot_resources/CID44971_scale.csv', row.names = 1)
data$Score = pmax(pmin(data[[var]], 1.75), -1.75)
data$Pathology = factor(data$Pathology,
                        levels = c('Normal + stroma + lymphocytes',
                                   'Stroma + adipose tissue',
                                   'Invasive cancer + lymphocytes',
                                   'Lymphocytes', 'DCIS', 'Stroma'))

data$size_cat <- cut(
  data$X.log.p_adj,
  breaks = c(-Inf, 2, 5, 20, 50, Inf),
  labels = c("ns", ">2", ">5", ">20", ">50"),
  include.lowest = TRUE
)

data$size_cat = as.character(data$size_cat)
data$size_cat[is.na(data$size_cat)] <- "ns"
data$size_cat = factor(data$size_cat, levels = c("ns", ">2", ">5", ">20", ">50"))

mat = data %>%
  group_by(Protein, Pathology) %>%
  summarise(meanScore = mean(Score, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Pathology, values_from = meanScore) %>%
  column_to_rownames('Protein')

# cluster
hc = hclust(dist(mat), method = 'complete')
protein_order = hc$labels[hc$order]
data$Protein = factor(data$Protein, levels = protein_order)

hc_col = hclust(dist(t(mat)), method = 'complete')
pathology_order = hc_col$labels[hc_col$order]
data$Pathology = factor(data$Pathology, levels = pathology_order)

score_rng   <- range(data$Score, na.rm = TRUE)
rng_width   <- diff(score_rng)
lower_lim   <- score_rng[1] - 0.05 * rng_width
upper_lim   <- score_rng[2] + 0.05 * rng_width

custom_colors <- rev(c(
  "#b35806", "#e08214", "#fdb863", "#fee0b6", "#f7f7f7",
  "#d8daeb", "#b2abd2", "#8073ac", "#542788"
))

p = ggplot(data, aes(x = Protein, y = Pathology)) +
  geom_point(aes(size = size_cat, color = Score)) +
  scale_color_gradientn(
    'Protein Mean\nDiff (scaled)',
    colors = custom_colors,
    limits = c(lower_lim, upper_lim)
  ) +
  scale_size_manual(
    name = '-log(p_adj)',
    values = c(
      'ns'   = 1,
      ">2"   = 2,
      ">5"  = 3,
      ">20" = 4,
      ">50"= 5
    )
  ) +
  theme_classic() +
  ggtitle('CID44971 Dotplot') +
  theme(
    panel.grid       = element_blank(),
    strip.placement  = 'outside',
    strip.background = element_rect(colour='grey', fill='white'),
    strip.text.y.left= element_text(face='bold', size=rel(RELSIZE), angle=0),
    plot.title       = element_text(hjust=0.5, size=rel(RELSIZE)),
    legend.title     = element_text(size=rel(RELSIZE)),
    legend.text      = element_text(size=rel(RELSIZE)),
    legend.key.height= unit(0.4, 'cm'),
    axis.title       = element_blank(),
    axis.text.x      = element_text(angle=90, vjust=0.5, hjust=1),
    text             = element_text(size=FONTSIZE)
  )

pdf('./Figure_4E.pdf', width=12, height=4)
print(p)
dev.off()