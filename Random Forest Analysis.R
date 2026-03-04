# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

# ===================== Random Forest Analysis for Core ASVs (Compatible with CCA/RDA, iMeta Publication Level) =====================
rm(list = ls())
options(stringsAsFactors = FALSE, warn = -1)

# 1. Load required packages (reuse with existing code to avoid duplicate installation)
required_packages <- c("tidyverse", "ggplot2", "randomForest", "ggpubr")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
options(warn = 0)
set.seed(123)  # Fixed seed for reproducible results (consistent with volcano plot code)

# 2. Read data (reuse existing files, no additional data required)
otu_raw <- read.csv("", row.names = 1, check.names = FALSE)  # Reuse abundance table
group_raw <- read.csv("", row.names = 1, check.names = FALSE)  # Reuse group table

# 3. Data preprocessing (compatible with CCA/RDA format, eliminate zero values)
otu_rf <- otu_raw + 1e-7  # Consistent with volcano plot preprocessing logic to avoid data inconsistency
otu_rf_t <- t(otu_rf) %>% as.data.frame()  # Transpose: samples as rows, ASVs as columns (compatible with random forest)
otu_rf_t$group <- factor(group_raw$group[match(rownames(otu_rf_t), rownames(group_raw))])

# 4. Train random forest model (focus on group-driven core ASVs, aligned with CCA/RDA)
rf_model <- randomForest(
  group ~ .,
  data = otu_rf_t,
  ntree = 500,  # Balance computational speed and result stability (common parameter in iMeta)
  importance = TRUE,  # Calculate ASV importance
  proximity = FALSE
)

# 5. Extract Top15 core ASVs (avoid plot overcrowding, linked with CCA/RDA results)
importance_df <- importance(rf_model) %>%
  as.data.frame() %>%
  rownames_to_column("ASV_ID") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice(1:15)  # Show only Top15 to focus on core ASVs precisely

# 6. Generate publication-level random forest plot (unified color/style with existing plots)
rf_plot <- ggplot(importance_df, aes(x = reorder(ASV_ID, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#FF7F0E", width = 0.7) +  # Coordinated with existing plot color scheme
  coord_flip() +
  labs(
    x = "Archaeal ASV ID",
    y = "Mean Decrease in Gini Index",
    title = "Top 15 Core ASVs Driving Community Separation (Random Forest)"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),  # Unified font with volcano/CCA plots
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(angle = 0, hjust = 1, face = "italic")
  )

# 7. Save plots (300dpi, PDF+PNG, consistent with existing plot format)
print(rf_plot)
ggsave("Archaeal_ASV_RandomForest_Core_Species.pdf", rf_plot, width = 9, height = 7, dpi = 300, device = cairo_pdf)
ggsave("Archaeal_ASV_RandomForest_Core_Species.png", rf_plot, width = 9, height = 7, dpi = 300)

# 8. Save importance results (supplement to differential analysis table)
write.csv(importance_df, "Archaeal_ASV_RandomForest_Core_Species_Results.csv", row.names = FALSE)
