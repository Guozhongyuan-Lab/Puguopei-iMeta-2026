# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

if (!require("vegan")) install.packages("vegan")

library(readxl)
library(vegan)
library(ggplot2)

data <- read_excel("")
group <- read_excel("")
data <- data.frame(data)
group <- data.frame(group)

rownames(data) <- data[, 1]
data <- data[, -1]
data <- t(data)

# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(data, method = "bray")

# Perform PCoA analysis
pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)

# Create coordinate dataframe
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points$Sample <- rownames(pcoa_points)

# Add group information
pcoa_points$Group <- subset(pcoa_points, grepl("CON", Sample), select = Sample, drop = TRUE)
pcoa_points$Group <- ifelse(is.na(pcoa_points$Group), "LP", "CON")

# Extract explained variance percentage
explained_var <- round(100 * pcoa_result$eig / sum(pcoa_result$eig), 2)

# Plot PCoA figure
ggplot(pcoa_points, aes(x = -PCoA1, y = PCoA2, color = Group)) +
  geom_point(size = 8, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.90, linetype = 1, alpha = 1, linewidth = 1) +
  scale_color_manual(values = c("CON" = "#B07FBF", "LP" = "#FFE287")) +
  labs(
    title = "",
    x = paste0("PCoA1 (", explained_var[1], "%)"),
    y = paste0("PCoA2 (", explained_var[2], "%)"),
    color = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )