
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

library(readxl)
# install.packages(c("tidyverse", "ggalluvial","patchwork","ggsci","devtools"))
library(tidyverse)
library(ggalluvial)
library(patchwork)
# devtools::install_github("sz-zyp/ggstyle")
library(ggstyle)

# Read pathway data
df <- read_excel("")

# Data preprocessing - split metabolites column
df_long <- df %>%
  separate_rows(Metabolites, sep = "/") %>%
  mutate(Pathway = factor(Pathway)) %>%
  mutate(Metabolites = factor(Metabolites))

# Prepare sankey plot data
df_sankey <- to_lodes_form(df_long %>% select(c("Metabolites", "Pathway")),
                           key = "x",
                           axes = c(1,2)) %>%
  mutate(flow_color = rep(df_long$Pathway, 2))

# Define color palette
d3_colors <- c(
  "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
  "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
  "#CAB2D6", "#FFFF99", "#8DD3C7", "#BEBADA", "#FB8072",
  "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9",
  "#BC80BD", "#CCEBC5", "#FFED6F", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
)

# Generate sankey plot
sankey_plot <- ggplot(data = df_sankey,
                      aes(x = x,
                          stratum = factor(stratum, levels = unique(stratum)),
                          alluvium = alluvium,
                          y = 1,
                          label = stratum,
                          fill = stratum
                      )) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw(base_family = "serif") +
  geom_flow(aes(fill = flow_color), alpha = 0.3, width = 0, knot.pos = 0.1) +
  geom_stratum(width = 0.05, color = "white") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5, family = "serif",
            hjust = 1, nudge_x = -0.04) +
  guides(fill = FALSE, color = FALSE) +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0, 0, 0, 3), units = "cm")
  ) +
  scale_x_discrete(expand = c(0.4, 0, 0, 0)) +
  scale_fill_manual(values = d3_colors)

print(sankey_plot)

# Prepare bubble plot data
bubble_df <- df %>%
  mutate(Pathway = factor(Pathway, levels = rev(df$Pathway))) %>%
  arrange(Pathway) %>%
  mutate(Pathway_num = cumsum(Hits) - Hits / 2)

# Generate bubble plot
bubble_plot <- ggplot(bubble_df, 
                      aes(x = -log10(PValue),
                          y = Pathway_num, 
                          color = PValue)) +
  geom_point(aes(size = `Enrichment Ratio`)) +
  scale_x_continuous(
    limits = c(0, 3.3),
    expand = c(0, 0),
    breaks = seq(0, 4, by = 1)
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, sum(bubble_df$Hits, na.rm = T))) +
  scale_color_gradient(low = "#FF6F48", high = "#4ECB90") +
  scale_radius(
    range = c(1, 5),
    name = "Enrichment Ratio"
  ) +
  guides(
    color = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  ) +
  theme_bw(base_family = "serif") +
  labs(size = "Enrichment Ratio", color = "P-value", y = "", x = expression(-log["10"]("P-value"))) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "inches"),
    panel.grid = element_blank()
  )

# Combine plots
combined_plot <- sankey_plot + bubble_plot + plot_layout(widths = c(2.5, 1))
print(combined_plot)

# Save plots
ggsave(combined_plot, file = "", width = 40, height = 10, units = 'cm', dpi = 1200)
ggsave(combined_plot, file = "", width = 40, height = 10, units = 'cm', dpi = 1200)
ggsave(combined_plot, file = "", width = 40, height = 10, units = 'cm', dpi = 1200)