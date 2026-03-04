# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

rm(list = ls())

packages <- c("psych", "dplyr", "tidyr", "tibble", "readxl", "pheatmap", "RColorBrewer", "stringi")
installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages], dependencies = TRUE, repos = "https://cloud.r-project.org/")
}
lapply(packages, library, character.only = TRUE)

clean_text <- function(text) {
  text <- stri_trans_general(text, "latin-ascii")
  text <- gsub("[^A-Za-z0-9 \\-\\_\\(\\)\\.]", " ", text)
  text <- gsub("\\s+", " ", text)
  text <- trimws(text)
  ifelse(text == "", "unknown", text)
}

simplify_name <- function(name, max_length = 12) {
  if (grepl("^Bacteria_", name)) {
    name <- gsub("^Bacteria_", "Bact_", name)
  } else if (grepl("^Archaea_", name)) {
    name <- gsub("^Archaea_", "Arch_", name)
  } else if (grepl("^Eukaryote_", name)) {
    name <- gsub("^Eukaryote_", "Euk_", name)
  } else if (grepl("^Fungi_", name)) {
    name <- gsub("^Fungi_", "Fung_", name)
  }
  if (nchar(name) > max_length) {
    name <- substr(name, 1, max_length)
  }
  return(name)
}

microbe_files <- list(
  Bacteria = "",
  Archaea = "",
  Eukaryote = "",
  Fungi = ""
)

missing_files <- c()
for (type in names(microbe_files)) {
  file_path <- file.path("", microbe_files[[type]])
  if (!file.exists(file_path)) {
    file_path <- file.path("", microbe_files[[type]])
    if (!file.exists(file_path)) {
      missing_files <- c(missing_files, microbe_files[[type]])
    }
  }
}
if (length(missing_files) > 0) {
  stop(paste("Missing microbe files: ", paste(missing_files, collapse = ", "), 
             "\nPlease confirm file location: ", ""))
}

process_microbe_data <- function(file_path, microbe_type) {
  if (!file.exists(file_path)) {
    file_path <- file.path("", basename(file_path))
  }
  
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    data <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
  } else if (grepl("\\.xlsx$", file_path, ignore.case = TRUE)) {
    data <- read_excel(file_path, na = c("", "NA"))
  } else {
    stop(paste("Unsupported file format: ", file_path))
  }
  
  colnames(data) <- sapply(colnames(data), clean_text)
  data_t <- as.data.frame(t(data))
  colnames(data_t) <- data_t[1, ]
  data_t <- data_t[-1, , drop = FALSE]
  
  rownames(data_t) <- toupper(gsub("[^A-Za-z0-9]", "", rownames(data_t)))
  
  data_t[] <- lapply(data_t, function(x) {
    as.numeric(as.character(gsub(",", ".", x)))
  })
  
  colnames(data_t) <- paste0(microbe_type, "_", make.unique(colnames(data_t), sep = "_"))
  
  data_t <- data_t[, colSums(is.na(data_t)) < nrow(data_t)]
  
  return(data_t)
}

microbe_list <- lapply(names(microbe_files), function(type) {
  process_microbe_data(microbe_files[[type]], type)
})
names(microbe_list) <- names(microbe_files)
merged_microbe <- bind_cols(microbe_list)

merged_microbe <- merged_microbe[, !duplicated(colnames(merged_microbe))]
merged_microbe <- merged_microbe[, colSums(is.na(merged_microbe)) < nrow(merged_microbe)]

cat("Merged microbe matrix: ", nrow(merged_microbe), " samples × ", ncol(merged_microbe), " ASVs\n")
if (ncol(merged_microbe) == 0) stop("No valid microbe data! Please check input files")

metab_file <- ""
metab_file_path <- file.path("", metab_file)

if (!file.exists(metab_file_path)) {
  metab_file_path <- file.path("", metab_file)
  if (!file.exists(metab_file_path)) {
    stop(paste("Missing metabolite file: ", metab_file, "\nPlease confirm file location: ", ""))
  }
}

metab_matrix <- read.csv(metab_file_path, row.names = 1, stringsAsFactors = FALSE, 
                         check.names = FALSE, na.strings = c("", "NA"))

colnames(metab_matrix) <- make.unique(sapply(colnames(metab_matrix), clean_text), sep = "_")

rownames(metab_matrix) <- toupper(gsub("[^A-Za-z0-9]", "", rownames(metab_matrix)))

metab_matrix <- metab_matrix[rowSums(is.na(metab_matrix)) < ncol(metab_matrix), ]
metab_matrix <- metab_matrix[, colSums(is.na(metab_matrix)) < nrow(metab_matrix)]

cat("Metabolite matrix: ", nrow(metab_matrix), " samples × ", ncol(metab_matrix), " metabolites\n")
if (ncol(metab_matrix) == 0) stop("No valid metabolite data! Please check input files")

common_samples <- intersect(rownames(merged_microbe), rownames(metab_matrix))
cat("Number of common samples matched: ", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No overlapping samples! Please check sample name format:\n",
       "Microbe sample names: ", paste(head(rownames(merged_microbe), 5), collapse = ", "), "\n",
       "Metabolite sample names: ", paste(head(rownames(metab_matrix), 5), collapse = ", "))
}

merged_microbe <- merged_microbe[common_samples, , drop = FALSE]
metab_matrix <- metab_matrix[common_samples, , drop = FALSE]

total_samples <- nrow(merged_microbe)

filtered_microbes <- merged_microbe %>%
  select(where(~ mean(., na.rm = TRUE) > 0.01)) %>%
  select(where(~ sum(. > 0, na.rm = TRUE) / total_samples > 0.5))

if (ncol(filtered_microbes) > 100) {
  top_microbes <- names(sort(colMeans(filtered_microbes, na.rm = TRUE), decreasing = TRUE)[1:100])
  filtered_microbes <- filtered_microbes[, top_microbes]
  cat("Number of microbe ASVs > 100, retaining top 100 by abundance\n")
}

cat("Filtered microbes: ", nrow(filtered_microbes), " samples × ", ncol(filtered_microbes), " ASVs\n")
if (ncol(filtered_microbes) == 0) stop("No microbes meet criteria! Suggest lowering abundance threshold to 0.005%")

calc_metab_stats <- function(col) {
  non_zero_vals <- col[col > 0 & !is.na(col)]
  if (length(non_zero_vals) == 0) return(c(coverage = 0, cv = 0, mean_ab = 0))
  coverage <- length(non_zero_vals) / total_samples
  cv <- sd(non_zero_vals) / mean(non_zero_vals)
  mean_ab <- mean(non_zero_vals)
  return(c(coverage = coverage, cv = cv, mean_ab = mean_ab))
}
metab_stats <- as.data.frame(t(apply(metab_matrix, 2, calc_metab_stats)))

filtered_metabs <- metab_matrix %>%
  select(which(
    metab_stats$coverage >= 0.3 &
      metab_stats$cv >= 0.2 &
      metab_stats$mean_ab >= 1000
  ))

if (ncol(filtered_metabs) > 50) {
  top_metabs <- names(sort(colMeans(filtered_metabs, na.rm = TRUE), decreasing = TRUE)[1:50])
  filtered_metabs <- filtered_metabs[, top_metabs]
  cat("Number of metabolites > 50, retaining top 50 by abundance\n")
}

cat("Filtered metabolites: ", nrow(filtered_metabs), " samples × ", ncol(filtered_metabs), " metabolites\n")
if (ncol(filtered_metabs) == 0) stop("No metabolites meet criteria! Suggest lowering abundance threshold to 500")

cat("Calculating Spearman correlations (FDR correction)...\n")
cor_result <- corr.test(
  x = filtered_microbes, 
  y = filtered_metabs, 
  method = "spearman",
  adjust = "fdr",
  ci = FALSE
)

cor_matrix <- cor_result$r    
p_matrix <- cor_result$p      

cor_matrix[is.na(cor_matrix)] <- 0
p_matrix[is.na(p_matrix)] <- 1

cat("Correlation calculation completed!\n")

cor_thresh <- 0.5    
p_thresh <- 0.1     

sig_indices <- which(abs(cor_matrix) >= cor_thresh & p_matrix <= p_thresh, arr.ind = TRUE)

if (nrow(sig_indices) == 0) {
  cat("No associations meet |r|≥0.5, P≤0.1, relaxing to |r|≥0.4, P≤0.15\n")
  cor_thresh <- 0.4
  p_thresh <- 0.15
  sig_indices <- which(abs(cor_matrix) >= cor_thresh & p_matrix <= p_thresh, arr.ind = TRUE)
}

sig_cor_df <- data.frame(
  Microbe = rownames(cor_matrix)[sig_indices[, "row"]],
  Metabolite = colnames(cor_matrix)[sig_indices[, "col"]],
  Correlation = cor_matrix[sig_indices],
  FDR_P_value = p_matrix[sig_indices]
) %>%
  distinct(Microbe, Metabolite, .keep_all = TRUE) %>%
  arrange(desc(abs(Correlation)))

top_n <- min(25, nrow(sig_cor_df))
top_cor_df <- sig_cor_df[1:top_n, ]

cat("Final selection of TOP", top_n, "core significant associations\n")
if (top_n == 0) {
  stop("No significant associations! Suggest further lowering thresholds:\n",
       "- Microbe abundance: 0.005%\n",
       "- Correlation: |r|≥0.3\n",
       "- P-value: ≤0.2")
}

top_microbes <- unique(top_cor_df$Microbe)
top_metabs <- unique(top_cor_df$Metabolite)
top_cor_matrix <- cor_matrix[top_microbes, top_metabs, drop = FALSE]

rownames(top_cor_matrix) <- sapply(rownames(top_cor_matrix), simplify_name)
colnames(top_cor_matrix) <- sapply(colnames(top_cor_matrix), simplify_name)

microbe_type <- gsub("_.+", "", top_microbes)
microbe_annotation <- data.frame(Microbe_Type = factor(microbe_type))
rownames(microbe_annotation) <- sapply(top_microbes, simplify_name)

custom_colors <- list(
  Microbe_Type = c(
    Bacteria = "#2E86AB",
    Archaea = "#E63946",
    Eukaryote = "#F1A208",
    Fungi = "#592E83"
  )[unique(microbe_type)]
)

cat("\nMicrobe type color mapping:\n")
for (type in names(custom_colors$Microbe_Type)) {
  cat(type, " → ", custom_colors$Microbe_Type[type], "\n")
}

heatmap_file <- file.path("", "TOP25_Microbe_Metabolite_Correlation_Heatmap.pdf")
heatmap_file <- normalizePath(heatmap_file, mustWork = FALSE)

pheatmap(
  top_cor_matrix,
  main = "TOP25 Significant Spearman Correlations (Microbes vs C/N Metabolites)",
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  breaks = seq(-1, 1, length.out = 100),
  annotation_row = microbe_annotation,
  annotation_colors = custom_colors,
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 8,
  treeheight_col = 8,
  fontsize = 6,
  fontsize_row = 5.5,
  fontsize_col = 5.5,
  angle_col = 45,
  display_numbers = FALSE,
  number_format = "%.2f",
  number_color = "black",
  filename = heatmap_file,
  width = 16,
  height = 12,
  dpi = 300,
  border_color = NA,
  legend = TRUE,
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  legend_labels = c("-1", "-0.5", "0", "0.5", "1")
)

cat("TOP25 heatmap saved to: ", heatmap_file, "\n")

top_cor_df <- top_cor_df %>%
  mutate(
    Microbe_Simplified = sapply(Microbe, simplify_name),
    Metabolite_Simplified = sapply(Metabolite, simplify_name),
    Microbe_Type = gsub("_.+", "", Microbe)
  ) %>%
  relocate(Microbe_Type, .after = Microbe)

result_files <- list(
  top_cor = "TOP25_Microbe_Metabolite_Correlations.csv",
  full_cor = "All_Significant_Correlations.csv",
  microbe_map = "Microbe_Name_Mapping.csv",
  metab_map = "Metabolite_Name_Mapping.csv",
  microbe_abund = "Filtered_Microbe_Abundance.csv",
  metab_abund = "Filtered_Metabolite_Abundance.csv"
)

write.csv(top_cor_df, file.path("", result_files$top_cor), row.names = FALSE)
write.csv(sig_cor_df, file.path("", result_files$full_cor), row.names = FALSE)
write.csv(
  data.frame(Microbe = top_microbes, Simplified_Name = sapply(top_microbes, simplify_name)),
  file.path("", result_files$microbe_map), row.names = FALSE
)
write.csv(
  data.frame(Metabolite = top_metabs, Simplified_Name = sapply(top_metabs, simplify_name)),
  file.path("", result_files$metab_map), row.names = FALSE
)
write.csv(filtered_microbes, file.path("", result_files$microbe_abund), row.names = TRUE)
write.csv(filtered_metabs, file.path("", result_files$metab_abund), row.names = TRUE)

ociation heatmap\n")