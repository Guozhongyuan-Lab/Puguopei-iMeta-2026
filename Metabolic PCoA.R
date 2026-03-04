# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

install.packages(c("ggplot2", "factoextra","ggrepel","openxlsx","vegan","patchwork"))
library(ggplot2)
library(ggrepel)
library(openxlsx)
library(vegan)
library(patchwork)

# Read metabolomics data
metabo_data <- read.xlsx("", sheet = 1, check.names = FALSE)

# Extract numerical matrix
value_matrix <- as.matrix(metabo_data[, 14:ncol(metabo_data)])
rownames(value_matrix) <- metabo_data$ID

# Create group information
sample_groups <- data.frame(
  sample = c("all_PSS1","all_PSS2","all_PSS3","all_PSS4","all_PSS5","all_PSS6",
             "all_NPSS1","all_NPSS2","all_NPSS3","all_NPSS4","all_NPSS5","all_NPSS6"),
  group = c(rep("PSS", 6), rep("NPSS", 6))
)
# Define plot parameters
plot_colors <- c("#2f7eb0", "#d00000")
point_shapes <- c(16, 16, 16)
legend_cols <- 1

# PCoA analysis
dist_matrix <- vegdist(t(value_matrix), method = 'bray')
pcoa <- cmdscale(dist_matrix, k = 5, eig = TRUE)
pc_points <- data.frame(pcoa$points[, 1:2])
colnames(pc_points) <- c('dim1', 'dim2')

# Calculate explained variance
pc_variance <- round(pcoa$eig/sum(pcoa$eig)*100, 2)[1:2]

# Merge group information
plot_data <- data.frame(
  sample = rownames(pc_points),
  dim1 = pc_points$dim1,
  dim2 = pc_points$dim2,
  group = sample_groups$group
)

# PERMANOVA analysis
adonis_result <- adonis2(dist_matrix ~ group, data = sample_groups, permutations = 999)
R2_value <- adonis_result$R2[1]
p_val <- adonis_result$`Pr(>F)`[1]
adonis_label <- sprintf("PERMANOVA:\nR² = %.3f\nP = %.3f", R2_value, p_val)

# Pairwise PERMANOVA
dist_mat <- as.matrix(dist_matrix)
group_levels <- unique(sample_groups$group)
combos <- combn(group_levels, 2, simplify = FALSE)

pairwise_results <- list()
for (i in seq_along(combos)) {
  g1 <- combos[[i]][1]
  g2 <- combos[[i]][2]
  
  idx <- sample_groups$group %in% c(g1, g2)
  sub_mat <- dist_mat[idx, idx]
  sub_dist <- as.dist(sub_mat)
  sub_group <- factor(sample_groups$group[idx])
  
  res <- adonis2(sub_dist ~ sub_group, permutations = 999)
  
  pairwise_results[[i]] <- data.frame(
    Comparison = paste(g1, "vs", g2),
    R2 = round(res$R2[1], 3),
    p_value = round(res$`Pr(>F)`[1], 3)
  )
}

pairwise_df <- do.call(rbind, pairwise_results)
print(pairwise_df)

# Custom theme
custom_theme <- function(){
  theme(
    text = element_text(family = "serif"),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 12),
    legend.position = "right"
  )
}

# Main plot
main_plot <- ggplot(plot_data, aes(dim1, dim2)) +
  geom_point(aes(color = group, shape = group), size = 3) +
  stat_ellipse(aes(color = group), level = 0.95) +
  scale_color_manual(values = plot_colors) +
  scale_shape_manual(values = point_shapes) +
  labs(
    x = sprintf("PCoA1 (%s%%)", pc_variance[1]),
    y = sprintf("PCoA2 (%s%%)", pc_variance[2])
  ) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    panel.grid = element_blank(),
    axis.title = element_text(color = 'black', family = "serif", size =15, margin = margin(t = 0, r = 0, b =-100, l = -100)),
    axis.text = element_text(colour = 'black', family = "serif",size =13, margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.ticks = element_line(color = 'black'),
    axis.ticks.length = unit(0, "lines"),
    legend.title = element_blank(),
    legend.text = element_text(family = "serif", size = 10),
    legend.key = element_blank(),
    legend.position = c(0.75, 0.80),
    legend.background = element_blank()
  ) +
  guides(col = guide_legend(ncol = legend_cols), shape = guide_legend(ncol = legend_cols))

print(main_plot)

# Marginal boxplot 1 (dim1)
boxplot_dim1 <- ggplot(plot_data, aes(x = group, y = dim1, fill = group)) +
  geom_boxplot(show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.35, size = 0.5) +
  geom_jitter(show.legend = FALSE) +
  scale_fill_manual(values = plot_colors) +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = "black", size = 15, family = "serif")
  )
print(boxplot_dim1)

# Marginal boxplot 2 (dim2)
boxplot_dim2 <- ggplot(plot_data, aes(x = group, y = dim2, fill = group)) +
  geom_boxplot(width = 0.5, show.legend = FALSE) +
  stat_boxplot(geom = "errorbar", width = 0.35, size = 0.5) +
  geom_jitter(width = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = plot_colors) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_line(color = "black"),
    axis.text.x = element_text(colour = "black", family = "serif", size = 15, angle = 45, vjust = 1, hjust = 1)
  )
print(boxplot_dim2)

# Statistical label plot
stat_label_plot <- ggplot(plot_data,aes(dim1,dim2)) +
  geom_text(aes(x=-0.5,y=0.6, label = adonis_label), size=4, family="serif") +
  theme_bw() +
  xlab("") + ylab("") +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
print(stat_label_plot)

# Plot layout 1 (deprecated, kept for compatibility)
temp_plot <- boxplot_dim1 + stat_label_plot + main_plot + boxplot_dim2 +
  plot_layout(heights=c(3,7),widths=c(7,3),ncol=2,nrow=2)
print(temp_plot)

# Final plot layout
final_plot <- (boxplot_dim1 + stat_label_plot) / (main_plot + boxplot_dim2) +
  plot_layout(heights = c(2, 5), widths = c(5, 2))
print(final_plot)

# Save plot
ggsave(final_plot, file = "", width = 13, height = 10, units = 'cm', dpi = 1200)
ggsave(final_plot, file = "", width = 13, height = 10, units = 'cm', dpi = 1200)
ggsave(final_plot, file = "", width = 13, height = 10, units = 'cm', dpi = 1200)