# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

# Script 2: C/N cycling related microbe-metabolite correlation analysis in PE plastisphere
# Focus: Key factors related to C/N cycling and plastisphere, not random selection

# Clear workspace
rm(list = ls())

# --- Step 0: Load required packages and configure environment ---
packages <- c("dplyr", "pheatmap", "showtext", "tibble", "stringr", "tidyr", "psych", "RColorBrewer")
installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages], dependencies = TRUE, repos = "https://cloud.r-project.org/")
}
lapply(packages, library, character.only = TRUE)

# Chinese font configuration (compatible with multiple systems)
library(showtext)
if (Sys.info()["sysname"] == "Windows") {
  font_add("SimHei", "C:/Windows/Fonts/simhei.ttf") 
} else if (Sys.info()["sysname"] == "Darwin") {
  font_add("SimHei", "/Library/Fonts/SimHei.ttf")
} else {
  font_add("SimHei", "/usr/share/fonts/truetype/liberation/simhei.ttf")
}
showtext_auto()
showtext_opts(dpi = 300) 

cat("--- Running PE plastisphere C/N cycling microbe-metabolite analysis script ---\n")

# --- Step 1: Data loading and ASV-to-Genus aggregation ---
# 1.1 Load ASV abundance data
abundance_file <- ""
if (!file.exists(abundance_file) && abundance_file != "") {
  stop(paste("❌ ASV abundance table not found: ", abundance_file, "\nPlease check file path"))
} else if (abundance_file != "") {
  asv_abundance <- read.csv(abundance_file, row.names = 1, check.names = FALSE)
} else {
  stop("❌ Please specify the ASV abundance file path")
}

# 1.2 Load taxonomy mapping table
tax_file <- ""
if (!file.exists(tax_file) && tax_file != "") {
  stop(paste("❌ Taxonomy mapping table not found: ", tax_file, "\nRun script 1 first to generate this file"))
} else if (tax_file != "") {
  taxonomy <- read.csv(tax_file, check.names = FALSE)
  taxonomy <- taxonomy %>% select(ASV_ID, Genus)
} else {
  stop("❌ Please specify the taxonomy file path")
}

# 1.3 Aggregate ASVs to genus level
genus_abundance <- asv_abundance %>%
  rownames_to_column("Sample") %>%
  tidyr::pivot_longer(cols = -Sample, names_to = "ASV_ID", values_to = "Abundance") %>%
  inner_join(taxonomy, by = "ASV_ID") %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus %in% c("", "unclassified", "Unclassified"), 
                        "Unclassified_genus", Genus)) %>%
  group_by(Sample, Genus) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE), .groups = 'drop') %>%
  tidyr::pivot_wider(names_from = Genus, values_from = Total_Abundance, values_fill = 0) %>%
  column_to_rownames("Sample")

cat(paste("✅ Data aggregation completed, obtained", ncol(genus_abundance), "genera\n"))

# 1.4 Load metabolite data
metab_file <- ""
if (metab_file != "") {
  if (!file.exists(metab_file)) {
    metab_file_alt <- ""
    if (!file.exists(metab_file_alt)) {
      stop(paste("❌ Metabolite table not found, provide one of the following files:\n",
                 "- Filtered_Metabolite_Abundance.csv\n",
                 "- cn_metab_matrix_final_clean.csv"))
    } else {
      metab_file <- metab_file_alt
    }
  }
  metab_abundance <- read.csv(metab_file, row.names = 1, check.names = FALSE)
} else {
  stop("❌ Please specify the metabolite file path")
}

# --- Step 2: Sample alignment ---
common_samples <- intersect(rownames(genus_abundance), rownames(metab_abundance))
if (length(common_samples) < 3) {
  stop(paste("❌ Only", length(common_samples), "common samples, cannot perform correlation analysis"))
}
genus_abundance <- genus_abundance[common_samples, , drop = FALSE]
metab_abundance <- metab_abundance[common_samples, , drop = FALSE]
cat(paste("✅ Sample alignment completed, valid samples: ", length(common_samples), "\n"))

# --- Step 3: Filter key factors based on research objectives (core optimization) ---
# 3.1 Filter core microbes related to C/N cycling + PE plastisphere
# Selection logic:
# - Nitrifying bacteria: Candidatus_Nitrocosmicus/Nitrososphaera (N cycling - nitrification)
# - Heterotrophic degraders: Halomonas (plastisphere dominant, C degradation), Escherichia-Shigella (N metabolism)
# - Other C/N cycling bacteria: Massilia (organic matter degradation)
target_genus_keywords <- c("Nitrocosmicus", "Nitrososphaera", "Halomonas", "Escherichia-Shigella", "Massilia")
key_genera_list <- colnames(genus_abundance)[sapply(colnames(genus_abundance), function(x) {
  any(sapply(target_genus_keywords, function(k) substr(x, 1, nchar(k)) == k))
})]

# 3.2 Filter core metabolites related to C/N cycling
# Selection logic:
# - N cycling: Uric acid, L Glutamic acid, L Glutamate, L Ornithine
# - C metabolism: Betaine, Choline, Leucine, Valine
# - C/N co-metabolism: Nicotinamide, Proline
target_metab_keywords <- c("Betaine", "Choline", "Uric acid", "L Glutamic acid", "L Glutamate", 
                           "Leucine", "Valine", "L Ornithine", "Nicotinamide", "Proline")
key_metab_list <- colnames(metab_abundance)[sapply(colnames(metab_abundance), function(x) {
  any(sapply(target_metab_keywords, function(k) grepl(k, x, ignore.case = TRUE)))
})]

# 3.3 Print filtering results (verify compliance with research objectives)
cat("\n===== Key factors filtered by research objectives =====\n")
cat("📌 C/N cycling related microbes (", length(key_genera_list), "):\n")
print(key_genera_list)
cat("\n📌 C/N cycling core metabolites (", length(key_metab_list), "):\n")
print(key_metab_list)

# 3.4 Non-empty validation (ensure valid factors are filtered)
if (length(key_genera_list) == 0) {
  stop("❌ No C/N cycling related microbes filtered, check keywords or data")
}
if (length(key_metab_list) == 0) {
  stop("❌ No C/N cycling related metabolites filtered, check keywords or data")
}

# 3.5 Trim data matrices
key_genera_data <- genus_abundance[, key_genera_list, drop = FALSE]
key_metab_data <- metab_abundance[, key_metab_list, drop = FALSE]

# --- Step 4: Calculate correlation matrix (Spearman + FDR correction for ecological data) ---
cor_test_result <- psych::corr.test(
  key_genera_data, 
  key_metab_data, 
  method = "spearman",  # Non-parametric test preferred for ecological data
  adjust = "fdr",       # FDR correction to control false positives (ecological standard)
  ci = FALSE           
)

correlation_matrix <- cor_test_result$r
p_matrix <- cor_test_result$p 

# Replace NA values (avoid plotting errors)
correlation_matrix[is.na(correlation_matrix)] <- 0
p_matrix[is.na(p_matrix)] <- 1

# --- Step 5: Plot scientific correlation heatmap (suitable for publication) ---
# Significance marking function (ecological standard)
star_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***",
                ifelse(p < 0.01, "**",
                       ifelse(p < 0.05, "*", "")
                )
         )
  )
}
p_labels <- matrix(star_p(p_matrix), nrow = nrow(p_matrix), ncol = ncol(p_matrix))

# Color-blind friendly palette (journal standard, avoid red-green)
color_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)

# Save heatmap (high resolution for publication)
heatmap_file <- "PE_plastisphere_CN_cycling_microbe_metabolite_correlation.pdf"
pdf(heatmap_file, width = 14, height = 10, family = "Helvetica") # Journal standard font
pheatmap(
  correlation_matrix,
  main = "Correlation between C/N Cycling Microbes and Metabolites in PE Plastisphere",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = p_labels,
  number_color = "black",
  number_fontsize = 9,
  fontsize_row = 11,
  fontsize_col = 11,
  color = color_palette,
  breaks = seq(-1, 1, length.out = 100),
  border_color = NA,
  cellwidth = 30, 
  cellheight = 25, 
  angle_col = 45,
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 15,
  treeheight_col = 15,
  legend = TRUE,
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  legend_labels = c("-1", "-0.5", "0", "0.5", "1"),
  annotation_legend = TRUE
)
dev.off()

# --- Step 6: Save analysis results (facilitate subsequent mechanism interpretation) ---
# Organize correlation results (include coefficient, P-value, significance)
cor_result_df <- as.data.frame(correlation_matrix) %>%
  rownames_to_column("Microbe_Genus") %>%
  pivot_longer(cols = -Microbe_Genus, names_to = "Metabolite", values_to = "Spearman_Correlation") %>%
  left_join(
    as.data.frame(p_matrix) %>%
      rownames_to_column("Microbe_Genus") %>%
      pivot_longer(cols = -Microbe_Genus, names_to = "Metabolite", values_to = "FDR_P_Value"),
    by = c("Microbe_Genus", "Metabolite")
  ) %>%
  mutate(Significance = star_p(FDR_P_Value)) %>%
  arrange(Microbe_Genus, desc(abs(Spearman_Correlation))) # Sort by absolute correlation value

# Save key factor list (with annotation of selection basis)
key_factors_df <- data.frame(
  Microbe = key_genera_list,
  Microbe_Function = c(
    "Nitrifying bacteria (N cycling - ammonia oxidation)",
    "Nitrifying bacteria (N cycling - ammonia oxidation)",
    "Plastisphere heterotrophic bacteria (C degradation + N metabolism)",
    "Enteric bacteria (N metabolism)",
    "Organic matter degrading bacteria (C cycling)"
  )[1:length(key_genera_list)],
  Metabolite = key_metab_list,
  Metabolite_Function = c(
    "Betaine (C/N co-metabolism)",
    "Choline (C/N metabolism)",
    "Leucine (C metabolism)",
    "L-Glutamic acid (N cycling)",
    "Proline (C/N metabolism)",
    "Uric acid (N cycling)",
    "Valine (C metabolism)",
    "Ornithine (N cycling)",
    "Nicotinamide (C/N metabolism)"
  )[1:length(key_metab_list)]
)

# Save all result files
write.csv(cor_result_df, "PE_plastisphere_CN_cycling_correlation_results.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(key_factors_df, "CN_cycling_key_factors_function_annotation.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(genus_abundance, "Genus_level_abundance_table_PE_plastisphere.csv", row.names = TRUE, fileEncoding = "UTF-8")
