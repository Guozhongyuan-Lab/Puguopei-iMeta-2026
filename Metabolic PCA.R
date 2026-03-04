# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

install.packages(c("ggplot2", "factoextra","ggrepel"))
install.packages("openxlsx")

library(ggplot2)
library(factoextra)
library(ggrepel)
library(openxlsx)

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

# Calculate row variance and filter
row_variance <- apply(value_matrix, 1, var, na.rm = TRUE)
filtered_matrix <- value_matrix[row_variance > 0, ]

# Perform PCA analysis
pca_result <- prcomp(t(filtered_matrix), scale. = TRUE)

# Extract PCA scores
pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$Group <- sample_groups$group

# Calculate explained variance
variance_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

# Plot PCA
PCA_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  stat_ellipse(
    type = "t",
    level = 0.95,
    geom = "path",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  scale_color_manual(values = c("#2f7eb0","#d00000")) +
  labs(
    x = paste0("PC1 (", variance_explained[1], "%)"),
    y = paste0("PC2 (", variance_explained[2], "%)"),
    color = "Group"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "serif"),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA),
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", colour = "black", size = 10),
    axis.text = element_text(face = "bold", colour = "black", size = 10),
    plot.title = element_text(face = "bold", colour = "black", size = 10, hjust = 0.5),
    legend.text = element_text(face = "bold", size = 10, color = "black"),
    legend.title = element_blank()
  )

# Preview plot
print(PCA_plot)

# Save plot
ggsave(PCA_plot, file = "", width = 13, height = 10, units = 'cm', dpi = 1200)
ggsave(PCA_plot, file = "", width = 13, height = 10, units = 'cm', dpi = 1200)
ggsave(PCA_plot, file = "", width = 13, height = 10, units = 'cm', dpi = 1200)