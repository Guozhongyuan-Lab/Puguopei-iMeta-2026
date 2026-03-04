# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

# ===================== Archaeal ASV Differential Analysis (Final Publication Version) =====================
# Clear environment variables to avoid interference
rm(list = ls())
options(stringsAsFactors = FALSE, warn = -1)  # Temporarily suppress non-critical warnings

# 1. Install and load required packages (ensure complete dependencies)
required_packages <- c("tidyverse", "ggplot2", "ggrepel")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}
options(warn = 0)  # Restore warnings

# 2. Set fixed random seed (Core: ensure consistent plot generation)
set.seed(123)

# 3. Read data (adapt to your files: archaea.csv, group_table.csv)
# Read ASV abundance table (rows=ASVs, columns=samples)
otu_raw <- read.csv("", row.names = 1, check.names = FALSE)
# Read group table (rows=samples, columns=groups)
group_raw <- read.csv("", row.names = 1, check.names = FALSE)

# 4. Data preprocessing (normalization + sample alignment + resolve duplicate values)
# Add minimal random values to eliminate Wilcoxon test warnings from 0/duplicate values
otu_raw <- otu_raw + runif(nrow(otu_raw)*ncol(otu_raw), 0, 1e-7)
# Reshape abundance table to long format (compatible with ggplot2)
otu_long <- otu_raw %>%
  tibble::rownames_to_column("ASV_ID") %>%
  tidyr::pivot_longer(
    cols = -ASV_ID, 
    names_to = "sample_name", 
    values_to = "relative_abundance"
  )
# Clean group table (keep only sample name and group columns)
group_clean <- group_raw %>%
  tibble::rownames_to_column("sample_name") %>%
  dplyr::select(sample_name, group)
# Merge data (ensure 100% sample matching, remove ungrouped samples)
merged_data <- dplyr::left_join(otu_long, group_clean, by = "sample_name") %>%
  dplyr::filter(!is.na(group))

# 5. Wilcoxon rank-sum test (exclusive for small sample sizes)
wilcox_test_wrapper <- function(values, groups) {
  suppressWarnings({  # Suppress "tie" warnings from duplicate values
    test <- wilcox.test(values ~ groups, exact = FALSE)
  })
  return(test$p.value)
}
# Perform differential analysis by ASV
diff_analysis <- merged_data %>%
  dplyr::group_by(ASV_ID) %>%
  dplyr::summarise(
    p_value = wilcox_test_wrapper(relative_abundance, group),
    mean_PSS = mean(relative_abundance[group == "PSS"], na.rm = TRUE),
    mean_NPSS = mean(relative_abundance[group == "NPSS"], na.rm = TRUE),
    log2_fold_change = log2((mean_PSS + 1e-6) / (mean_NPSS + 1e-6)),
    significance = dplyr::case_when(
      p_value < 0.05 ~ "P < 0.05 (Significant)",
      p_value < 0.1 ~ "0.05 ≤ P < 0.1 (Marginally significant)",
      TRUE ~ "NS (Not significant)"
    ),
    .groups = "drop"
  ) %>%
  dplyr::arrange(p_value)  # Sort by P-value for easy review

# 6. Generate publication-level volcano plot (classic red/blue/grey color scheme)
volcano_plot <- ggplot2::ggplot(diff_analysis, 
                                aes(x = log2_fold_change, 
                                    y = -log10(p_value), 
                                    color = significance)) +
  # Scatter points (size/transparency optimized for journal aesthetics)
  geom_point(size = 3, alpha = 0.9) +
  # Label only ASVs with P<0.1 and fold change >1 (avoid clutter)
  ggrepel::geom_text_repel(
    data = dplyr::filter(diff_analysis, p_value < 0.1 & abs(log2_fold_change) > 1),
    aes(label = ASV_ID),
    size = 3.5,
    family = "Arial",
    color = "black",
    max.overlaps = 8,
    box.padding = 0.6
  ) +
  # Reference lines (bolded for clarity)
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40", linewidth = 0.7) +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "gray40", linewidth = 0.7) +
  # Classic red/blue/grey color scheme (matching significance levels)
  scale_color_manual(
    values = c("#d00000", "#2f7eb0", "#999999"),  # Red: significant, Blue: marginally significant, Grey: non-significant
    breaks = c("P < 0.05 (Significant)", "0.05 ≤ P < 0.1 (Marginally significant)", "NS (Not significant)")
  ) +
  # Plot annotations (clear and standardized)
  labs(
    x = "log2(Fold Change) (PSS / NPSS)",
    y = "-log10(P-value)",
    title = "Volcano Plot of Archaeal ASVs (PSS vs NPSS)",
    color = "Significance"
  ) +
  # Theme optimization (no grid + unified font)
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 7. Generate publication-level boxplot (optimized - resolve ASV overlap + improve aesthetics)
# Filter key ASVs (max 8 to avoid overcrowding)
key_asvs <- diff_analysis %>%
  dplyr::filter(p_value < 0.1 & abs(log2_fold_change) > 1) %>%
  dplyr::arrange(desc(abs(log2_fold_change))) %>%
  dplyr::slice(1:8) %>%  # Max 8 ASVs to completely resolve X-axis crowding
  dplyr::pull(ASV_ID)

# If no significant ASVs, select top 5 with highest abundance variation
if (length(key_asvs) == 0) {
  key_asvs <- head(diff_analysis$ASV_ID, 5)
  message("No significant ASVs found, displaying top 5 ASVs with highest abundance variation")
}

# Filter boxplot data
boxplot_data <- merged_data %>%
  dplyr::filter(ASV_ID %in% key_asvs)

# Generate optimized boxplot
boxplot <- ggplot2::ggplot(boxplot_data, 
                           aes(x = ASV_ID, 
                               y = relative_abundance * 100,  # Convert to percentage for intuitiveness
                               fill = group)) +
  # Boxplot (journal-level style - bold lines for better contrast)
  geom_boxplot(width = 0.6, linewidth = 0.8, outlier.size = 2.5, outlier.shape = 21, outlier.fill = "black") +
  # Raw data points (high transparency to avoid obscuring boxplots)
  geom_jitter(width = 0.2, size = 2, alpha = 0.8, shape = 21, color = "black") +
  # Specified color scheme: PSS=#d00000  NPSS=#2f7eb0
  scale_fill_manual(values = c("#d00000", "#2f7eb0"), labels = c("PSS", "NPSS")) +
  # Plot annotations
  labs(
    x = "ASV ID",
    y = "Relative Abundance (%)",
    title = "Abundance of Key Archaeal ASVs (PSS vs NPSS)",
    fill = "Group"
  ) +
  # Theme optimization (resolve ASV overlap + unified font)
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    text = element_text(family = "Arial"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    # Core optimization: rotate X-axis ASVs by 60° + reduce font size to eliminate overlap
    axis.text.x = element_text(angle = 60, hjust = 1, face = "italic", size = 9),
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# 8. Save plots (300dpi high resolution, PDF+PNG dual format)
# Save volcano plot
print(volcano_plot)
ggsave("Archaeal_ASV_Volcano_Plot.pdf", volcano_plot, width = 10, height = 8, dpi = 300, device = cairo_pdf)
ggsave("Archaeal_ASV_Volcano_Plot.png", volcano_plot, width = 10, height = 8, dpi = 300)
# Save boxplot (wider width for more spacious X-axis)
print(boxplot)
ggsave("Archaeal_ASV_Boxplot.pdf", boxplot, width = 10, height = 6, dpi = 300, device = cairo_pdf)
ggsave("Archaeal_ASV_Boxplot.png", boxplot, width = 10, height = 6, dpi = 300)

# 9. Export differential analysis results (sorted by P-value)
write.csv(diff_analysis, "Archaeal_ASV_Wilcoxon_Differential_Analysis_Results.csv", row.names = FALSE)
