# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

rm(list=ls())

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)
library(caret)
library(ggrepel)

# Read data
archaea <- read.csv("", row.names = 1, check.names = FALSE)
bacteria <- read.csv("", row.names = 1, check.names = FALSE)
eukaryote <- read.csv("", row.names = 1, check.names = FALSE)
fungi <- read.csv("", row.names = 1, check.names = FALSE)
soil <- read.csv("", row.names = 1, check.names = FALSE)

# Function: Perform RDA analysis and plotting
perform_rda <- function(community_data, soil_data, title, genus_column = NULL, 
                        label_distance = 0.5, label_size = 3, max_overlaps = 10) {
  community_t <- t(community_data)
  
  # Check and handle multicollinearity
  cor_matrix <- cor(soil_data)
  high_cor <- findCorrelation(cor_matrix, cutoff = 0.85)
  
  if (length(high_cor) > 0) {
    cat("Warning: Highly correlated variables detected, removing:\n")
    cat(names(soil_data)[high_cor], "\n")
    soil_data <- soil_data[, -high_cor]
  }
  
  # RDA analysis
  rda_model <- rda(community_t ~ ., data = soil_data)
  
  # Variance decomposition
  variance_explained <- RsquareAdj(rda_model)$adj.r.squared * 100
  
  # Extract RDA results
  rda_scores <- scores(rda_model, display = c("sites", "species", "bp"))
  
  # Create RDA plot
  p <- ggplot() +
    geom_point(data = as.data.frame(rda_scores$sites), 
               aes(x = RDA1, y = RDA2), 
               size = 3, color = "#2f7eb0", alpha = 0.7) +
    geom_segment(data = as.data.frame(rda_scores$biplot), 
                 aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2), 
                 arrow = arrow(length = unit(0.2, "cm")), 
                 color = "#E09D94") +
    geom_text(data = as.data.frame(rda_scores$biplot), 
              aes(x = RDA1 * 2.2, y = RDA2 * 2.2, label = rownames(rda_scores$biplot)), 
              color = "#E09D94", size = 4)
  
  # Add species labels
  if (!is.null(genus_column)) {
    top_species <- rownames(rda_scores$species[1:10, ])
    top_genera <- genus_column[match(top_species, rownames(community_data))]
    
    species_data <- as.data.frame(rda_scores$species[1:10, ])
    species_data$Genus <- top_genera
    
    p <- p +
      geom_point(data = species_data, 
                 aes(x = RDA1, y = RDA2), 
                 size = 3, color = "#FFB55F", shape = 17) +
      geom_text_repel(data = species_data, 
                      aes(x = RDA1, y = RDA2, label = Genus),
                      color = "#FFB55F", size = label_size,
                      box.padding = unit(label_distance, "lines"),
                      max.overlaps = max_overlaps,
                      segment.color = 'grey50',
                      segment.size = 0.3)
  } else {
    p <- p +
      geom_point(data = as.data.frame(rda_scores$species[1:10, ]), 
                 aes(x = RDA1, y = RDA2), 
                 size = 3, color = "#FFB55F", shape = 17) +
      geom_text_repel(data = as.data.frame(rda_scores$species[1:10, ]), 
                      aes(x = RDA1, y = RDA2, label = rownames(rda_scores$species[1:10, ])),
                      color = "#FFB55F", size = label_size,
                      box.padding = unit(label_distance, "lines"),
                      max.overlaps = max_overlaps,
                      segment.color = 'grey50',
                      segment.size = 0.3)
  }
  
  # Set plot labels
  p <- p +
    labs(title = paste0(title, " RDA Analysis (Adj R²: ", 
                        round(variance_explained, 2), "%)"),
         x = paste0("RDA1 (", round(rda_model$CCA$eig[1]/sum(rda_model$CCA$eig)*100, 2), "%)"),
         y = paste0("RDA2 (", round(rda_model$CCA$eig[2]/sum(rda_model$CCA$eig)*100, 2), "%)")) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          legend.position = "bottom")
  
  return(list(model = rda_model, plot = p))
}

# Export RDA results to CSV
export_rda_results_csv <- function(rda_result, soil_data, title, dir_prefix, genus_data = NULL) {
  dir.create(dir_prefix, showWarnings = FALSE)
  
  env_contribution <- scores(rda_result$model, display = "bp")
  write.csv(env_contribution, file = file.path(dir_prefix, "env_contribution.csv"))
  
  species_scores <- scores(rda_result$model, display = "species")
  
  if (!is.null(genus_data)) {
    species_with_genus <- data.frame(
      ASV = rownames(species_scores),
      Genus = genus_data[rownames(species_scores)],
      species_scores
    )
    write.csv(species_with_genus, file = file.path(dir_prefix, "species_scores.csv"))
  } else {
    write.csv(species_scores, file = file.path(dir_prefix, "species_scores.csv"))
  }
  
  site_scores <- scores(rda_result$model, display = "sites")
  write.csv(site_scores, file = file.path(dir_prefix, "site_scores.csv"))
  
  eigenvalues <- rda_result$model$CCA$eig
  variance_percent <- eigenvalues / sum(eigenvalues) * 100
  cum_variance_percent <- cumsum(variance_percent)
  
  variance_table <- data.frame(
    Axis = paste0("RDA", 1:length(eigenvalues)),
    Eigenvalue = eigenvalues,
    Variance = variance_percent,
    Cumulative = cum_variance_percent
  )
  
  write.csv(variance_table, file = file.path(dir_prefix, "variance.csv"))
  
  rda_summary <- summary(rda_result$model)
  rda_coef <- rda_summary$coefficients
  write.csv(as.data.frame(rda_coef), file = file.path(dir_prefix, "coefficients.csv"))
  
  cat(paste0("Successfully exported ", title, " RDA results to directory: ", dir_prefix, "\n"))
}

# Process microbial data for RDA
process_microbial_data <- function(microbial_data, soil_data, title) {
  genus <- microbial_data[, ncol(microbial_data)]
  names(genus) <- rownames(microbial_data)
  
  data <- microbial_data[, -ncol(microbial_data)]
  
  common_samples <- intersect(colnames(data), rownames(soil_data))
  data_filtered <- data[, common_samples]
  soil_filtered <- soil_data[common_samples, ]
  
  rda_result <- perform_rda(data_filtered, soil_filtered, title, genus, 
                            label_distance = 0.8, label_size = 3.5, max_overlaps = 20)
  
  export_rda_results_csv(rda_result, soil_filtered, title, paste0(title, "_rda_results"), genus)
  
  return(rda_result)
}

# Re-read data (removed duplicate read.csv calls)
archaea <- read.csv("", row.names = 1, check.names = FALSE)
bacteria <- read.csv("", row.names = 1, check.names = FALSE)
eukaryote <- read.csv("", row.names = 1, check.names = FALSE)
fungi <- read.csv("", row.names = 1, check.names = FALSE)
soil <- read.csv("", row.names = 1, check.names = FALSE)

# Process each microbial group
archaea_rda <- process_microbial_data(archaea, soil, "Archaea")
bacteria_rda <- process_microbial_data(bacteria, soil, "Bacteria")
eukaryote_rda <- process_microbial_data(eukaryote, soil, "Eukaryote")
fungi_rda <- process_microbial_data(fungi, soil, "Fungi")

# Print individual RDA plots
print(archaea_rda$plot)
print(bacteria_rda$plot)
print(eukaryote_rda$plot)
print(fungi_rda$plot)

# Merge microbial data
archaea_data <- archaea[, -ncol(archaea)]
archaea_genus <- archaea[, ncol(archaea)]
bacteria_data <- bacteria[, -ncol(bacteria)]
bacteria_genus <- bacteria[, ncol(bacteria)]
eukaryote_data <- eukaryote[, -ncol(eukaryote)]
eukaryote_genus <- eukaryote[, ncol(eukaryote)]
fungi_data <- fungi[, -ncol(fungi)]
fungi_genus <- fungi[, ncol(fungi)]

all_microbes <- rbind(
  data.frame(ASV = rownames(archaea_data), Domain = "Archaea", archaea_data, Genus = archaea_genus),
  data.frame(ASV = rownames(bacteria_data), Domain = "Bacteria", bacteria_data, Genus = bacteria_genus),
  data.frame(ASV = rownames(eukaryote_data), Domain = "Eukaryote", eukaryote_data, Genus = eukaryote_genus),
  data.frame(ASV = rownames(fungi_data), Domain = "Fungi", fungi_data, Genus = fungi_genus)
)

# Calculate domain abundance
domain_abundance <- all_microbes %>%
  group_by(Domain) %>%
  summarise(across(where(is.numeric), sum)) %>%
  tibble::column_to_rownames(var = "Domain")

domain_abundance_t <- t(domain_abundance)
domain_abundance_t <- domain_abundance_t[rownames(soil), ]

# Check sample and variable counts
n_samples <- nrow(soil)
n_vars <- ncol(soil)

cat("Number of samples:", n_samples, "\n")
cat("Original number of variables:", n_vars, "\n")

# Strict variable filtering
cor_matrix <- cor(soil)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.8)

if (length(high_cor) > 0) {
  cat("Warning: Highly correlated variables detected, removing:\n")
  cat(names(soil)[high_cor], "\n")
  soil_filtered <- soil[, -high_cor]
} else {
  soil_filtered <- soil
}

# Remove near zero variance variables
nzv <- nearZeroVar(soil_filtered)
if (length(nzv) > 0) {
  cat("Warning: Near zero variance variables detected, removing:\n")
  cat(names(soil_filtered)[nzv], "\n")
  soil_filtered <- soil_filtered[, -nzv]
}

# Ensure variable count < sample count
max_vars <- max(5, floor(n_samples/3))
filtered_vars <- ncol(soil_filtered)

if (filtered_vars > max_vars) {
  cat("Warning: Too many variables, further reduction...\n")
  if (require(randomForest, quietly = TRUE)) {
    rf_model <- randomForest(x = soil_filtered, y = as.factor(rownames(domain_abundance_t)), ntree = 500)
    var_importance <- importance(rf_model)
    top_vars <- order(var_importance[,1], decreasing = TRUE)[1:max_vars]
    soil_for_rda <- soil_filtered[, top_vars]
    cat("Variables selected by random forest for RDA:\n")
    cat(names(soil_for_rda), "\n")
  } else {
    variances <- apply(soil_filtered, 2, var)
    top_vars <- order(variances, decreasing = TRUE)[1:max_vars]
    soil_for_rda <- soil_filtered[, top_vars]
    cat("Variables selected by variance for RDA:\n")
    cat(names(soil_for_rda), "\n")
  }
} else {
  soil_for_rda <- soil_filtered
  cat("Variables used for RDA:\n")
  cat(names(soil_for_rda), "\n")
}

# Final check and RDA analysis
final_vars <- ncol(soil_for_rda)
if (final_vars >= n_samples - 1) {
  cat("Error: Variable count close to/exceeds sample count.\n")
  cat("Consider collecting more samples or reducing variables.\n")
} else {
  cat("Starting RDA analysis...\n")
  all_rda <- rda(domain_abundance_t ~ ., data = soil_for_rda)
  
  if (!is.null(all_rda$CCA) && !is.null(attr(all_rda$CCA, "rank"))) {
    if (attr(all_rda$CCA, "rank") == 0) {
      cat("Error: RDA model overfitted, cannot proceed.\n")
    } else {
      all_variance_explained <- RsquareAdj(all_rda)$adj.r.squared * 100
      all_scores <- scores(all_rda, display = c("sites", "species", "bp"))
      
      if (is.null(all_scores$biplot) || nrow(all_scores$biplot) == 0) {
        cat("Warning: No biplot data for environmental factors.\n")
        empty_biplot <- data.frame(RDA1 = numeric(0), RDA2 = numeric(0))
        
        all_plot <- ggplot() +
          geom_point(data = as.data.frame(all_scores$sites), 
                     aes(x = RDA1, y = RDA2), 
                     size = 3, color = "#93C3C6", alpha = 0.7) +
          geom_point(data = as.data.frame(all_scores$species), 
                     aes(x = RDA1, y = RDA2, color = rownames(all_scores$species)), 
                     size = 5, shape = 16) +
          geom_text(data = as.data.frame(all_scores$species), 
                    aes(x = RDA1, y = RDA2, label = rownames(all_scores$species)), 
                    color = "black", size = 4, vjust = -1.2) +
          labs(title = paste0("All Microbes RDA Analysis (Adj R²: ", 
                              round(all_variance_explained, 2), "%)"),
               x = paste0("RDA1 (", round(all_rda$CCA$eig[1]/sum(all_rda$CCA$eig)*100, 2), "%)"),
               y = paste0("RDA2 (", round(all_rda$CCA$eig[2]/sum(all_rda$CCA$eig)*100, 2), "%)"),
               color = "Microbial Domain") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.position = "bottom",
                legend.title = element_text(face = "bold"))
      } else {
        species_scores <- as.data.frame(all_scores$species)
        species_scores$Domain <- rownames(species_scores)
        
        all_plot <- ggplot() +
          geom_point(data = as.data.frame(all_scores$sites), 
                     aes(x = RDA1, y = RDA2), 
                     size = 3, color = "#93C3C6", alpha = 0.7) +
          geom_segment(data = as.data.frame(all_scores$biplot), 
                       aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2), 
                       arrow = arrow(length = unit(0.2, "cm")), 
                       color = "#E09D94") +
          geom_text(data = as.data.frame(all_scores$biplot), 
                    aes(x = RDA1 * 2.2, y = RDA2 * 2.2, label = rownames(all_scores$biplot)), 
                    color = "#E09D94", size = 4) +
          geom_point(data = species_scores, 
                     aes(x = RDA1, y = RDA2, color = Domain), 
                     size = 5, shape = 16) +
          geom_text(data = species_scores, 
                    aes(x = RDA1, y = RDA2, label = Domain), 
                    color = "black", size = 4, vjust = -1.2) +
          labs(title = paste0("All Microbes RDA Analysis (Adj R²: ", 
                              round(all_variance_explained, 2), "%)"),
               x = paste0("RDA1 (", round(all_rda$CCA$eig[1]/sum(all_rda$CCA$eig)*100, 2), "%)"),
               y = paste0("RDA2 (", round(all_rda$CCA$eig[2]/sum(all_rda$CCA$eig)*100, 2), "%)"),
               color = "Microbial Domain") +
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.position = "bottom",
                legend.title = element_text(face = "bold"))
      }
      
      print(all_plot)
      export_rda_results_csv(list(model = all_rda), soil_for_rda, "All Microbes", "all_microbes_rda_results")
    }
  } else {
    cat("Error: RDA model failed to fit properly.\n")
  }
}