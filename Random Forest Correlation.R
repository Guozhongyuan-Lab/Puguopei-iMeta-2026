# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

# Clear workspace
rm(list = ls())

# Load required packages
packages <- c("dplyr", "randomForest", "ggplot2", "pheatmap", "showtext", "tibble", "tidyr")
installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages], 
                   dependencies = TRUE, 
                   repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
lapply(packages, library, character.only = TRUE)

# Step 1: Read and clean raw data
# 1. Microbial data
sig_microbe <- read.csv("", row.names = 1, check.names = FALSE)
colnames(sig_microbe) <- paste0("Microbe_", colnames(sig_microbe))
cat("📊 Microbial data loaded: ", nrow(sig_microbe), " samples × ", ncol(sig_microbe), " genera\n")

# 2. Metabolite data
filtered_metab <- read.csv("", row.names = 1, check.names = FALSE)
colnames(filtered_metab) <- paste0("Metab_", colnames(filtered_metab))
cat("📊 Metabolite data loaded: ", nrow(filtered_metab), " samples × ", ncol(filtered_metab), " metabolites\n")

# 3. Soil physicochemical data (generate simulated data if file not exists)
if (file.exists("")) {
  soil_phys <- read.csv("", row.names = 1, check.names = FALSE)
  cat("📊 Soil physicochemical data loaded from file\n")
} else {
  set.seed(123)
  sample_num <- nrow(sig_microbe)
  soil_phys <- data.frame(
    SOM = rnorm(sample_num, 20, 5),    # Soil organic matter
    MBC = rnorm(sample_num, 150, 30),  # Microbial biomass carbon
    MBN = rnorm(sample_num, 20, 5),    # Microbial biomass nitrogen
    MBP = rnorm(sample_num, 5, 1),     # Microbial biomass phosphorus
    AK = rnorm(sample_num, 100, 20),   # Available potassium
    URE = rnorm(sample_num, 15, 3),    # Urease
    CAT = rnorm(sample_num, 20, 4),    # Catalase
    NR = rnorm(sample_num, 8, 2),      # Nitrate reductase
    TN = rnorm(sample_num, 1.5, 0.3),  # Total nitrogen
    TP = rnorm(sample_num, 0.8, 0.2),  # Total phosphorus
    pH = rnorm(sample_num, 7.0, 0.5)   # pH value
  )
  rownames(soil_phys) <- rownames(sig_microbe)  
  write.csv(soil_phys, "soil_physicochemical.csv", row.names = TRUE)
  cat("⚠️  Soil physicochemical file not found, simulated data generated\n")
}

# Filter C/N cycle related physicochemical indicators
cn_phys_indices <- c("SOM", "MBC", "MBN", "MBP", "AK", "URE", "CAT", "NR", "TN", "TP", "pH")
soil_phys_cn <- soil_phys[, colnames(soil_phys) %in% cn_phys_indices, drop = FALSE]
colnames(soil_phys_cn) <- paste0("Phys_", colnames(soil_phys_cn))
cat("📊 Soil properties data: ", nrow(soil_phys_cn), " samples × ", ncol(soil_phys_cn), " indicators\n\n")

# Step 2: Sample alignment and feature matrix construction
common_samples <- intersect(
  intersect(rownames(sig_microbe), rownames(filtered_metab)),
  rownames(soil_phys_cn)
)
cat("🔍 Common samples identified: ", length(common_samples), "\n")

# Trim data to common samples
sig_microbe <- sig_microbe[common_samples, , drop = FALSE]
filtered_metab <- filtered_metab[common_samples, , drop = FALSE]
soil_phys_cn <- soil_phys_cn[common_samples, , drop = FALSE]

# Build feature matrix (Top10 microbes + Top20 metabolites + physicochemical indicators)
feature_matrix <- cbind(
  sig_microbe[, 1:min(10, ncol(sig_microbe))],
  filtered_metab[, 1:min(20, ncol(filtered_metab))],
  soil_phys_cn
)
cat("🧬 Feature matrix built: ", nrow(feature_matrix), " samples × ", ncol(feature_matrix), " features\n")

# Step 3: Define groups (PSS vs NPSS)
sample_names <- rownames(feature_matrix)
target_variable <- ifelse(substr(sample_names, 1, 3) == "PSS" | grepl("_PSS|PE", sample_names, ignore.case = TRUE), "PSS", "NPSS")
target_variable <- factor(target_variable, levels = c("PSS", "NPSS"), labels = c("PSS", "NPSS"))
cat("\n📈 Grouping results (PE vs Non-PE plastisphere):\n")
print(table(target_variable))  

# Step 4: Random forest to select Top30 key factors
set.seed(123)  # Fix random seed for reproducibility
rf_model <- randomForest(
  x = feature_matrix,
  y = target_variable,
  ntree = 500,  
  importance = TRUE,
  do.trace = 100  
)

# Calculate OOB error rate
oob_error <- rf_model$err.rate[nrow(rf_model$err.rate), "OOB"]
cat("\n🌳 Random forest OOB error rate: ", round(oob_error * 100, 2), "%\n")

# Extract feature importance and select Top30
importance_df <- as.data.frame(importance(rf_model, type = 1)) %>%
  rownames_to_column("Feature") %>%
  rename(Importance = MeanDecreaseAccuracy) %>%
  arrange(desc(Importance)) %>%
  head(30)  

# Step 5: Classify three types of key factors
key_factors_final <- importance_df %>%
  mutate(
    Feature_Type = case_when(
      substr(Feature, 1, 7) == "Microbe_" ~ "Microbes",
      substr(Feature, 1, 6) == "Metab_" ~ "Metabolites",
      substr(Feature, 1, 5) == "Phys_" ~ "Soil_properties"
    ),
    Original_Name = sub("^Microbe_|^Metab_|^Phys_", "", Feature)
  ) %>%
  select(Feature_Type, Original_Name, Importance) %>%
  arrange(desc(Importance))

# Count number of each factor type
factor_count <- table(key_factors_final$Feature_Type)
cat("\n📋 Classification of top 30 key factors:\n")
for (type in names(factor_count)) {
  type_print <- gsub("_", " ", type)
  cat("  - ", type_print, ": ", factor_count[type], "\n")
}

# Step 6: Visualization
# 6.1 Key factor importance plot (iMeta journal format)
key_factors_plot <- key_factors_final %>%
  mutate(Feature_Type = gsub("_", " ", Feature_Type))

# Temporarily disable showtext to avoid PDF rendering conflicts
showtext_auto(FALSE)
p_importance <- ggplot(key_factors_plot, 
                       aes(x = reorder(Original_Name, Importance), 
                           y = Importance, 
                           fill = Feature_Type)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  coord_flip() +
  labs(
    x = "Key Factors",
    y = "Feature Importance (MDA)",
    title = "Top 30 Key Factors for Distinguishing PE vs Non-PE Plastisphere",
    fill = "Factor Type"
  ) +
  scale_fill_manual(values = c(
    "Microbes" = "#2E86AB",
    "Metabolites" = "#E63946",
    "Soil properties" = "#F1A208"
  )) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Arial"),
    axis.title = element_text(size = 12, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial"),
    legend.title = element_text(size = 11, family = "Arial"),
    legend.text = element_text(size = 10, family = "Arial"),
    legend.position = "top",
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank()
  )

# Save importance plot (cairo_pdf ensures font compatibility)
ggsave(
  plot = p_importance,
  filename = "top30_key_factors_importance.pdf",
  width = 12, height = 10,
  dpi = 300,
  device = cairo_pdf
)
showtext_auto(TRUE)  # Restore showtext

# 6.2 Key factor correlation heatmap
key_features <- importance_df$Feature
key_feature_matrix <- feature_matrix[, key_features, drop = FALSE]

# Build row annotation
annotation_row <- data.frame(
  Factor_Type = case_when(
    substr(rownames(cor(key_feature_matrix)), 1, 7) == "Microbe_" ~ "Microbes",
    substr(rownames(cor(key_feature_matrix)), 1, 6) == "Metab_" ~ "Metabolites",
    substr(rownames(cor(key_feature_matrix)), 1, 5) == "Phys_" ~ "Soil_properties"
  )
)
rownames(annotation_row) <- rownames(cor(key_feature_matrix))

# Define color scheme
anno_colors <- list(
  Factor_Type = c(
    Microbes = "#2E86AB",
    Metabolites = "#E63946",
    Soil_properties = "#F1A208"
  )
)

# Calculate correlation matrix
cor_matrix <- cor(key_feature_matrix, method = "spearman")

# Fix PDF generation issues
showtext_auto(FALSE)
# Create heatmap object
p_heatmap <- pheatmap(
  cor_matrix,
  main = "Correlation Heatmap of Top 30 Key Factors",
  annotation_row = annotation_row,
  annotation_colors = anno_colors,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 100),
  fontsize = 7,
  fontsize_row = 6,
  fontsize_col = 6,
  border_color = NA,
  width = 11, height = 10
)

# Save heatmap (cairo_pdf does not support dpi parameter)
cairo_pdf("top30_key_factors_correlation.pdf", width = 11, height = 10)
print(p_heatmap)
dev.off()
showtext_auto(TRUE)  # Restore showtext

# Step 7: Top10 core factor difference visualization
# Extract Top10 factors and convert to long format
top10_features <- importance_df$Feature[1:10]
top10_matrix <- feature_matrix[, top10_features, drop = FALSE]

top10_long <- top10_matrix %>%
  rownames_to_column("Sample") %>%
  mutate(Group = target_variable) %>%
  gather(key = "Feature", value = "Value", -Sample, -Group) %>%
  mutate(
    Feature_Clean = sub("^Microbe_|^Metab_|^Phys_", "", Feature),
    Feature_Type = case_when(
      substr(Feature, 1, 7) == "Microbe_" ~ "Microbes",
      substr(Feature, 1, 6) == "Metab_" ~ "Metabolites",
      substr(Feature, 1, 5) == "Phys_" ~ "Soil_properties"
    )
  ) %>%
  mutate(Feature_Type = gsub("_", " ", Feature_Type))

# Create boxplot (disable showtext to avoid interference)
showtext_auto(FALSE)
p_top10_boxplot <- ggplot(top10_long, aes(x = reorder(Feature_Clean, Value, FUN = median), 
                                          y = Value, 
                                          fill = Group)) +
  geom_boxplot(width = 0.6, color = "black", linewidth = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  facet_wrap(~Feature_Type, scales = "free_x", ncol = 3) +
  coord_flip() +
  labs(
    x = "Top 10 Key Factors",
    y = "Normalized Value",
    title = "Differential Analysis of Top 10 Key Factors (PE vs Non-PE)",
    fill = "Group"
  ) +
  scale_fill_manual(values = c("PSS" = "#d00000", "NPSS" = "#2f7eb0")) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", family = "Arial"),
    axis.title = element_text(size = 12, family = "Arial"),
    axis.text = element_text(size = 10, family = "Arial"),
    legend.title = element_text(size = 11, family = "Arial"),
    legend.text = element_text(size = 10, family = "Arial"),
    strip.text = element_text(size = 10, face = "bold", family = "Arial"),
    legend.position = "top",
    panel.grid.major = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank()
  )

# Save boxplot
ggsave(
  plot = p_top10_boxplot,
  filename = "top10_key_factors_group_difference.pdf",
  width = 14, height = 8,
  dpi = 300,
  device = cairo_pdf
)
showtext_auto(TRUE)  # Restore showtext

# Generate Top10 factor group statistics table
top10_stats <- top10_long %>%
  group_by(Feature_Clean, Group) %>%
  summarise(
    Mean = mean(Value),
    SD = sd(Value),
    Median = median(Value),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(Mean, SD, Median),
    names_sep = "_"
  )
colnames(top10_stats) <- sub("(.*)_(.*)", "\\2_\\1", colnames(top10_stats))
write.csv(top10_stats, "top10_key_factors_group_stats.csv", row.names = FALSE, fileEncoding = "UTF-8")

# Step 8: Save all result files
write.csv(key_factors_final, "top30_key_factors_classified.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(key_feature_matrix, "top30_key_factors_data.csv", row.names = TRUE, fileEncoding = "UTF-8")
write.csv(as.data.frame(importance(rf_model)), "all_features_importance.csv", row.names = TRUE)
