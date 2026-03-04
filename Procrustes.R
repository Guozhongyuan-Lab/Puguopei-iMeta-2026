# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

rm(list=ls())

options(repos = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
install.packages(c("vegan", "ggplot2", "dplyr", "ggrepel", "caret"))

library(vegan)    
library(ggplot2)  
library(dplyr)    
library(ggrepel)  
library(caret)    

# Define iMeta plotting theme
theme_imeta <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 10),
      axis.text = element_text(size = 9),
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8),
      plot.margin = margin(10, 10, 10, 10, "mm")
    )
}

# Read microbial data function
read_microbial_data <- function(file_path, min_samples = 3) {
  cat(paste0("Reading microbial data: ", file_path, "\n"))
  df <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  numeric_cols <- sapply(df, is.numeric)
  if (sum(numeric_cols) == 0) {
    stop("Error: No numeric columns in microbial data, check data format")
  }
  df <- df[, numeric_cols, drop = FALSE]
  cat(paste0("Retained ", sum(numeric_cols), " numeric columns (removed ", sum(!numeric_cols), " non-numeric columns)\n"))
  
  present <- rowSums(df > 0) >= min_samples
  if (sum(!present) > 0) {
    cat(paste0("Filtered ", sum(!present), " low abundance species\n"))
    df <- df[present, , drop = FALSE]
  }
  
  df_t <- as.data.frame(t(df))
  rownames(df_t) <- make.names(rownames(df_t), unique = TRUE)
  cat(paste0("Successfully loaded ", nrow(df_t), " samples, ", ncol(df_t), " species\n"))
  
  if (any(is.na(df_t))) {
    cat("Warning: Data contains missing values, replacing with 0\n")
    df_t[is.na(df_t)] <- 0
  }
  
  return(df_t)
}

# Read soil data function
read_soil_data <- function(file_path) {
  cat(paste0("Reading soil data: ", file_path, "\n"))
  df <- read.csv(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  numeric_cols <- sapply(df, is.numeric)
  df <- df[, numeric_cols, drop = FALSE]
  rownames(df) <- make.names(rownames(df), unique = TRUE)
  cat(paste0("Successfully loaded ", nrow(df), " samples, ", ncol(df), " soil properties\n"))
  
  if (any(is.na(df))) {
    cat("Warning: Soil data contains missing values, replacing with column means\n")
    df <- as.data.frame(lapply(df, function(x) {
      ifelse(is.na(x), mean(x, na.rm = TRUE), x)
    }))
  }
  
  return(df)
}

# Procrustes analysis main function
perform_procrustes <- function(community_data, soil_data, title, 
                               method = "bray", transform = "hellinger", 
                               permutations = 999, save_plot = TRUE) {
  if (!is.numeric(community_data) || !is.matrix(community_data)) {
    community_data <- as.matrix(community_data)
  }
  if (!is.numeric(soil_data) || !is.matrix(soil_data)) {
    soil_data <- as.matrix(soil_data)
  }
  
  common_samples <- intersect(rownames(community_data), rownames(soil_data))
  if (length(common_samples) < 2) {
    stop("Error: Insufficient common samples (<2), cannot perform analysis")
  }
  cat(paste0("Found ", length(common_samples), " common samples, starting analysis...\n"))
  
  community_filtered <- community_data[common_samples, , drop = FALSE]
  soil_filtered <- soil_data[common_samples, , drop = FALSE]
  
  if (ncol(soil_filtered) > 1) {
    cor_matrix <- cor(soil_filtered)
    high_cor <- caret::findCorrelation(cor_matrix, cutoff = 0.9)
    if (length(high_cor) > 0) {
      cat("Warning: Removing highly correlated soil properties: ", paste(colnames(soil_filtered)[high_cor], collapse = ", "), "\n")
      soil_filtered <- soil_filtered[, -high_cor, drop = FALSE]
    }
  }
  
  cat("Checking community data...\n")
  if (any(is.na(community_filtered))) {
    cat("Warning: Community data contains missing values, replacing with 0\n")
    community_filtered[is.na(community_filtered)] <- 0
  }
  
  zero_var <- apply(community_filtered, 2, var) == 0
  if (sum(zero_var) > 0) {
    cat("Warning: Removing ", sum(zero_var), " zero variance species\n")
    community_filtered <- community_filtered[, !zero_var, drop = FALSE]
  }
  
  if (ncol(community_filtered) < 2) {
    stop("Error: Insufficient valid species (<2), cannot calculate distance matrix")
  }
  
  cat(paste0("Applying ", transform, " transformation...\n"))
  if (transform == "hellinger") {
    community_transformed <- decostand(community_filtered, method = "hellinger")
  } else if (transform == "log") {
    community_transformed <- log(community_filtered + 1)
  } else {
    community_transformed <- community_filtered
  }
  
  cat("Calculating community distance matrix...\n")
  community_dist <- tryCatch({
    vegdist(community_transformed, method = "bray")
  }, error = function(e) {
    stop("Distance calculation failed: ", e$message)
  })
  
  cat("Calculating soil data distance matrix...\n")
  soil_dist <- dist(soil_filtered, method = "euclidean")
  
  cat("Performing principal coordinate analysis...\n")
  k <- min(nrow(community_filtered) - 1, 2)  
  community_pcoa <- tryCatch({
    cmdscale(community_dist, k = k, eig = TRUE)
  }, error = function(e) {
    stop("Community data PCoA failed: ", e$message)
  })
  
  soil_pcoa <- tryCatch({
    cmdscale(soil_dist, k = k, eig = TRUE)
  }, error = function(e) {
    stop("Soil data PCoA failed: ", e$message)
  })
  
  if (is.null(community_pcoa$points) || ncol(community_pcoa$points) < 2) {
    stop("Error: Insufficient dimensions in community PCoA results")
  }
  if (is.null(soil_pcoa$points) || ncol(soil_pcoa$points) < 2) {
    stop("Error: Insufficient dimensions in soil PCoA results")
  }
  
  cat("Performing Procrustes analysis...\n")
  procrustes_result <- tryCatch({
    procrustes(community_pcoa$points, soil_pcoa$points)
  }, error = function(e) {
    stop("Procrustes analysis failed: ", e$message)
  })
  
  cat("Performing permutation test...\n")
  perm_test <- tryCatch({
    protest(community_pcoa$points, soil_pcoa$points, permutations = permutations)
  }, error = function(e) {
    stop("Permutation test failed: ", e$message)
  })
  
  rotated_community <- as.data.frame(community_pcoa$points)
  rotated_soil <- as.data.frame(soil_pcoa$points)
  colnames(rotated_community) <- c("PC1", "PC2")
  colnames(rotated_soil) <- c("PC1", "PC2")
  
  # Group identification (fixed logic)
  rotated_community$Sample <- rownames(rotated_community)
  rotated_community$Group <- case_when(
    substr(rotated_community$Sample, 1, 3) == "PSS" ~ "PSS",
    substr(rotated_community$Sample, 1, 4) == "NPSS" ~ "NPSS",
    TRUE ~ "Unknown"
  )
  
  rotated_soil$Sample <- rownames(rotated_soil)
  rotated_soil$Group <- case_when(
    substr(rotated_soil$Sample, 1, 3) == "PSS" ~ "PSS",
    substr(rotated_soil$Sample, 1, 4) == "NPSS" ~ "NPSS",
    TRUE ~ "Unknown"
  )
  
  # Print group identification results
  cat("Community data group identification results:\n")
  print(table(rotated_community$Group))
  cat("Soil data group identification results:\n")
  print(table(rotated_soil$Group))
  
  community_eig_perc <- round(community_pcoa$eig[1:2]/sum(community_pcoa$eig)*100, 1)
  soil_eig_perc <- round(soil_pcoa$eig[1:2]/sum(soil_pcoa$eig)*100, 1)
  
  # Define group colors
  group_colors <- c("PSS" = "#d00000", "NPSS" = "#2f7eb0", "Unknown" = "gray50")
  
  # Create plot
  p <- ggplot() +
    geom_point(data = rotated_community, aes(PC1, PC2, color = Group), 
               size = 3, shape = 16, alpha = 0.8, na.rm = TRUE) +
    geom_point(data = rotated_soil, aes(PC1, PC2, color = Group), 
               size = 3, shape = 17, alpha = 0.8, na.rm = TRUE) +
    geom_segment(data = data.frame(
      x = rotated_community$PC1, y = rotated_community$PC2,
      xend = rotated_soil$PC1, yend = rotated_soil$PC2,
      Sample = rotated_community$Sample,
      Group = rotated_community$Group
    ), aes(x, y, xend = xend, yend = yend, group = Sample), 
    linetype = "dashed", color = "gray50", alpha = 0.6, linewidth = 0.5,
    na.rm = TRUE) +
    geom_text_repel(data = rotated_community, 
                    aes(PC1, PC2, label = Sample, color = Group),
                    size = 3, max.overlaps = 50, box.padding = 0.5,
                    na.rm = TRUE) +
    geom_text_repel(data = rotated_soil, 
                    aes(PC1, PC2, label = Sample, color = Group),
                    size = 3, max.overlaps = 50, box.padding = 0.5,
                    na.rm = TRUE) +
    labs(title = title,
         subtitle = paste0("M² = ", round(procrustes_result$ss, 3), 
                           ", p = ", round(perm_test$signif, 3)),
         x = paste0("PC1 (", community_eig_perc[1], "%)"),
         y = paste0("PC2 (", community_eig_perc[2], "%)"),
         color = "Group") +
    scale_color_manual(values = group_colors) +
    theme_imeta() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  
  # Save plot
  if (save_plot) {
    png_filename <- paste0(title, "_procrustes_plot.png")
    ggsave(png_filename, plot = p, width = 180, height = 150, units = "mm", dpi = 300)
    tiff_filename <- paste0(title, "_procrustes_plot.tif")
    ggsave(tiff_filename, plot = p, width = 180, height = 150, units = "mm", dpi = 300, 
           device = "tiff", compression = "lzw")
    cat(paste0("Plot saved as: ", png_filename, " and ", tiff_filename, "\n"))
  }
  
  return(list(plot = p, result = procrustes_result, perm_test = perm_test,
              community_pcoa = community_pcoa, soil_pcoa = soil_pcoa))
}

# Data processing main workflow
process_data <- function(microbial_file, soil_data, title, 
                         method = "bray", transform = "hellinger") {
  cat(paste0("\n==== Processing ", title, " data ====\n"))
  microbial_data <- read_microbial_data(microbial_file, min_samples = 3)
  
  common <- intersect(rownames(microbial_data), rownames(soil_data))
  if (length(common) == 0) {
    stop("Error: No matching samples between microbial and soil data")
  }
  cat(paste0("Sample matching successful: ", length(common), " common samples\n"))
  
  proc_result <- tryCatch({
    perform_procrustes(microbial_data[common, , drop = FALSE], 
                       soil_data[common, , drop = FALSE], 
                       title, method, transform)
  }, error = function(e) {
    cat(paste0("Analysis failed: ", e$message, "\n"))
    return(NULL)
  })
  
  if (!is.null(proc_result)) {
    print(proc_result$plot)
    cat(paste0(title, " Procrustes analysis completed (M² = ", round(proc_result$result$ss, 3), 
               ", p = ", round(proc_result$perm_test$signif, 3), ")\n"))
  }
  
  return(proc_result)
}

# Execute analysis
soil <- read_soil_data("")
archaea_proc <- process_data("", soil, "Archaea")
bacteria_proc <- process_data("", soil, "Bacteria")
eukaryota_proc <- process_data("", soil, "Eukaryota")
fungi_proc <- process_data("", soil, "Fungi")

save(archaea_proc, bacteria_proc, eukaryota_proc, fungi_proc, 
     file = "all_procrustes_results.RData")
cat("All group analyses completed, results saved!\n")