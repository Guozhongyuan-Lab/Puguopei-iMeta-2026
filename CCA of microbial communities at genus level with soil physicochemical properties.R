# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

# Clear workspace
rm(list=ls())

# Load required packages
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(gridExtra)
library(caret)
library(ggrepel)

# Function: Perform CCA analysis and generate plot
perform_cca <- function(community_data, soil_data, title, genus_column = NULL, 
                        label_distance = 0.5, label_size = 3, max_overlaps = 10) {
  # Transpose community data (samples as rows, species as columns)
  community_t <- t(community_data)
  
  # Check and handle multicollinearity
  cor_matrix <- cor(soil_data)
  high_cor <- findCorrelation(cor_matrix, cutoff = 0.85)
  
  if (length(high_cor) > 0) {
    cat("Warning: Highly correlated variables detected, removing:\n")
    cat(names(soil_data)[high_cor], "\n")
    soil_data <- soil_data[, -high_cor]
  }
  
  # Perform CCA analysis
  cca_model <- cca(community_t ~ ., data = soil_data)
  
  # Variance decomposition
  variance_explained <- RsquareAdj(cca_model)$adj.r.squared * 100
  
  # Extract CCA results
  cca_scores <- scores(cca_model, display = c("sites", "species", "bp"))
  
  # Create CCA plot
  p <- ggplot() +
    # Add sample points
    geom_point(data = as.data.frame(cca_scores$sites), 
               aes(x = CCA1, y = CCA2), 
               size = 3, color = "#2f7eb0", alpha = 0.7) +
    # Add environmental factor arrows
    geom_segment(data = as.data.frame(cca_scores$biplot), 
                 aes(x = 0, y = 0, xend = CCA1 * 2, yend = CCA2 * 2), 
                 arrow = arrow(length = unit(0.2, "cm")), 
                 color = "#E09D94") +
    # Add environmental factor labels
    geom_text(data = as.data.frame(cca_scores$biplot), 
              aes(x = CCA1 * 2.2, y = CCA2 * 2.2, label = rownames(cca_scores$biplot)), 
              color = "#E09D94", size = 4)
  
  # Add species labels (Genus if provided, otherwise ASV)
  if (!is.null(genus_column)) {
    # Get top 10 species and their corresponding Genus
    top_species <- rownames(cca_scores$species[1:10, ])
    top_genera <- genus_column[match(top_species, rownames(community_data))]
    
    # Prepare species data with Genus labels
    species_data <- as.data.frame(cca_scores$species[1:10, ])
    species_data$Genus <- top_genera
    
    # Add species points with smart labels
    p <- p +
      geom_point(data = species_data, 
                 aes(x = CCA1, y = CCA2), 
                 size = 3, color = "#FFB55F", shape = 17) +
      geom_text_repel(data = species_data, 
                      aes(x = CCA1, y = CCA2, label = Genus),
                      color = "#FFB55F", size = label_size,
                      box.padding = unit(label_distance, "lines"),
                      max.overlaps = max_overlaps,
                      segment.color = 'grey50',
                      segment.size = 0.3)
  } else {
    # Add species points with ASV labels
    p <- p +
      geom_point(data = as.data.frame(cca_scores$species[1:10, ]), 
                 aes(x = CCA1, y = CCA2), 
                 size = 3, color = "#FFB55F", shape = 17) +
      geom_text_repel(data = as.data.frame(cca_scores$species[1:10, ]), 
                      aes(x = CCA1, y = CCA2, label = rownames(cca_scores$species[1:10, ])),
                      color = "#FFB55F", size = label_size,
                      box.padding = unit(label_distance, "lines"),
                      max.overlaps = max_overlaps,
                      segment.color = 'grey50',
                      segment.size = 0.3)
  }
  
  # Set plot titles and axis labels
  p <- p +
    labs(title = paste0(title, " CCA Analysis (Adj R²: ", 
                        round(variance_explained, 2), "%)"),
         x = paste0("CCA1 (", round(cca_model$CCA$eig[1]/sum(cca_model$CCA$eig)*100, 2), "%)"),
         y = paste0("CCA2 (", round(cca_model$CCA$eig[2]/sum(cca_model$CCA$eig)*100, 2), "%)")) +
    # Add theme settings
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          legend.position = "bottom")
  
  # Return CCA model and plot
  return(list(model = cca_model, plot = p))
}

# Function: Export CCA results to CSV files
export_cca_results_csv <- function(cca_result, soil_data, title, dir_prefix, genus_data = NULL) {
  # Create output directory
  dir.create(dir_prefix, showWarnings = FALSE)
  
  # 1. Export environmental factor contributions
  env_contribution <- scores(cca_result$model, display = "bp")
  write.csv(env_contribution, file = file.path(dir_prefix, "env_contribution.csv"))
  
  # 2. Export species scores (with Genus if available)
  species_scores <- scores(cca_result$model, display = "species")
  
  if (!is.null(genus_data)) {
    # Match species with Genus
    species_with_genus <- data.frame(
      ASV = rownames(species_scores),
      Genus = genus_data[rownames(species_scores)],
      species_scores
    )
    write.csv(species_with_genus, file = file.path(dir_prefix, "species_scores.csv"))
  } else {
    write.csv(species_scores, file = file.path(dir_prefix, "species_scores.csv"))
  }
  
  # 3. Export sample scores
  site_scores <- scores(cca_result$model, display = "sites")
  write.csv(site_scores, file = file.path(dir_prefix, "site_scores.csv"))
  
  # 4. Export eigenvalues and explained variance
  eigenvalues <- cca_result$model$CCA$eig
  variance_percent <- eigenvalues / sum(eigenvalues) * 100
  cum_variance_percent <- cumsum(variance_percent)
  
  variance_table <- data.frame(
    Axis = paste0("CCA", 1:length(eigenvalues)),
    Eigenvalue = eigenvalues,
    Variance = variance_percent,
    Cumulative = cum_variance_percent
  )
  
  write.csv(variance_table, file = file.path(dir_prefix, "variance.csv"))
  
  # 5. Export CCA model coefficients
  cca_summary <- summary(cca_result$model)
  cca_coef <- cca_summary$coefficients
  write.csv(as.data.frame(cca_coef), file = file.path(dir_prefix, "coefficients.csv"))
  
  cat(paste0("Successfully exported ", title, " CCA results to directory: ", dir_prefix, "\n"))
}

# Function: Process microbial data and perform CCA
process_microbial_data <- function(microbial_data, soil_data, title) {
  # Extract Genus column (last column)
  genus <- microbial_data[, ncol(microbial_data)]
  names(genus) <- rownames(microbial_data)
  
  # Remove Genus column (keep only sample data)
  data <- microbial_data[, -ncol(microbial_data)]
  
  # Filter common samples
  common_samples <- intersect(colnames(data), rownames(soil_data))
  data_filtered <- data[, common_samples]
  soil_filtered <- soil_data[common_samples, ]
  
  # Perform CCA analysis with adjusted label parameters
  cca_result <- perform_cca(data_filtered, soil_filtered, title, genus, 
                            label_distance = 0.8, label_size = 3.5, max_overlaps = 20)
  
  # Export results
  export_cca_results_csv(cca_result, soil_filtered, title, paste0(title, "_cca_results"), genus)
  
  return(cca_result)
}

# Read data
archaea <- read.csv("", row.names = 1, check.names = FALSE)
bacteria <- read.csv("", row.names = 1, check.names = FALSE)
eukaryote <- read.csv("", row.names = 1, check.names = FALSE)
fungi <- read.csv("", row.names = 1, check.names = FALSE)
soil <- read.csv("", row.names = 1, check.names = FALSE)

# Process each microbial group
archaea_cca <- process_microbial_data(archaea, soil, "Archaea")
bacteria_cca <- process_microbial_data(bacteria, soil, "Bacteria")
eukaryote_cca <- process_microbial_data(eukaryote, soil, "Eukaryote")
fungi_cca <- process_microbial_data(fungi, soil, "Fungi")

# Display individual CCA plots
print(archaea_cca$plot)
print(bacteria_cca$plot)
print(eukaryote_cca$plot)
print(fungi_cca$plot)

# Combine all microbial data
# Extract data and Genus information
archaea_data <- archaea[, -ncol(archaea)]
archaea_genus <- archaea[, ncol(archaea)]
bacteria_data <- bacteria[, -ncol(bacteria)]
bacteria_genus <- bacteria[, ncol(bacteria)]
eukaryote_data <- eukaryote[, -ncol(eukaryote)]
eukaryote_genus <- eukaryote[, ncol(eukaryote)]
fungi_data <- fungi[, -ncol(fungi)]
fungi_genus <- fungi[, ncol(fungi)]

# Merge microbial datasets
all_microbes <- rbind(
  data.frame(ASV = rownames(archaea_data), Domain = "Archaea", archaea_data, Genus = archaea_genus),
  data.frame(ASV = rownames(bacteria_data), Domain = "Bacteria", bacteria_data, Genus = bacteria_genus),
  data.frame(ASV = rownames(eukaryote_data), Domain = "Eukaryote", eukaryote_data, Genus = eukaryote_genus),
  data.frame(ASV = rownames(fungi_data), Domain = "Fungi", fungi_data, Genus = fungi_genus)
)

# Calculate domain-level relative abundance
domain_abundance <- all_microbes %>%
  group_by(Domain) %>%
  summarise(across(where(is.numeric), sum)) %>%
  tibble::column_to_rownames(var = "Domain")

# Transpose data (samples as rows, microbial domains as columns)
domain_abundance_t <- t(domain_abundance)

# Align sample order with soil data
domain_abundance_t <- domain_abundance_t[rownames(soil), ]

# Check sample and variable counts
n_samples <- nrow(soil)
n_vars <- ncol(soil)

cat("Number of samples:", n_samples, "\n")
cat("Original number of variables:", n_vars, "\n")

# Strict variable selection strategy
# 1. Handle multicollinearity in soil data
cor_matrix <- cor(soil)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.8)

if (length(high_cor) > 0) {
  cat("Warning: Highly correlated variables detected, removing:\n")
  cat(names(soil)[high_cor], "\n")
  soil_filtered <- soil[, -high_cor]
} else {
  soil_filtered <- soil
}

# 2. Remove near-zero variance variables
nzv <- nearZeroVar(soil_filtered)
if (length(nzv) > 0) {
  cat("Warning: Near-zero variance variables detected, removing:\n")
  cat(names(soil_filtered)[nzv], "\n")
  soil_filtered <- soil_filtered[, -nzv]
}

# 3. Ensure variable count is less than sample count
max_vars <- max(5, floor(n_samples/3))
filtered_vars <- ncol(soil_filtered)

if (filtered_vars > max_vars) {
  cat("Warning: Too many variables, further reduction needed...\n")
  # Variable selection using random forest or variance
  if (require(randomForest, quietly = TRUE)) {
    rf_model <- randomForest(x = soil_filtered, y = as.factor(rownames(domain_abundance_t)), ntree = 500)
    var_importance <- importance(rf_model)
    top_vars <- order(var_importance[,1], decreasing = TRUE)[1:max_vars]
    soil_for_cca <- soil_filtered[, top_vars]
    cat("Variables selected by random forest for CCA:\n")
    cat(names(soil_for_cca), "\n")
  } else {
    # Fallback: variable selection by variance
    variances <- apply(soil_filtered, 2, var)
    top_vars <- order(variances, decreasing = TRUE)[1:max_vars]
    soil_for_cca <- soil_filtered[, top_vars]
    cat("Variables selected by variance for CCA:\n")
    cat(names(soil_for_cca), "\n")
  }
} else {
  soil_for_cca <- soil_filtered
  cat("Variables used for CCA after filtering:\n")
  cat(names(soil_for_cca), "\n")
}

# Final check for valid CCA model
final_vars <- ncol(soil_for_cca)
if (final_vars >= n_samples - 1) {
  cat("Error: Variable count is close to or exceeds sample count.\n")
  cat("Please collect more samples or reduce variables further.\n")
} else {
  # Perform combined CCA analysis
  cat("Starting combined CCA analysis...\n")
  all_cca <- cca(domain_abundance_t ~ ., data = soil_for_cca)
  
  # Check if model is valid
  if (!is.null(all_cca$CCA) && !is.null(attr(all_cca$CCA, "rank"))) {
    if (attr(all_cca$CCA, "rank") == 0) {
      cat("Error: CCA model overfitted, cannot proceed with analysis.\n")
    } else {
      all_variance_explained <- RsquareAdj(all_cca)$adj.r.squared * 100
      all_scores <- scores(all_cca, display = c("sites", "species", "bp"))
      
      # Check if biplot data exists
      if (is.null(all_scores$biplot) || nrow(all_scores$biplot) == 0) {
        cat("Warning: Environmental factor biplot data is empty, cannot draw arrows.\n")
        # Create empty biplot dataframe
        empty_biplot <- data.frame(CCA1 = numeric(0), CCA2 = numeric(0))
        
        # Create combined CCA plot (without environmental factors)
        all_plot <- ggplot() +
          # Add sample points
          geom_point(data = as.data.frame(all_scores$sites), 
                     aes(x = CCA1, y = CCA2), 
                     size = 3, color = "#93C3C6", alpha = 0.7) +
          
          # Add microbial domain points (colored by domain)
          geom_point(data = as.data.frame(all_scores$species), 
                     aes(x = CCA1, y = CCA2, color = rownames(all_scores$species)), 
                     size = 5, shape = 16) +
          
          # Add microbial domain labels
          geom_text(data = as.data.frame(all_scores$species), 
                    aes(x = CCA1, y = CCA2, label = rownames(all_scores$species)), 
                    color = "black", size = 4, vjust = -1.2) +
          
          # Set plot titles and labels
          labs(title = paste0("All Microbes CCA Analysis (Adj R²: ", 
                              round(all_variance_explained, 2), "%)"),
               x = paste0("CCA1 (", round(all_cca$CCA$eig[1]/sum(all_cca$CCA$eig)*100, 2), "%)"),
               y = paste0("CCA2 (", round(all_cca$CCA$eig[2]/sum(all_cca$CCA$eig)*100, 2), "%)"),
               color = "Microbial Domain") +
          
          # Theme settings
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.position = "bottom",
                legend.title = element_text(face = "bold"))
      } else {
        # Prepare species scores with domain labels
        species_scores <- as.data.frame(all_scores$species)
        species_scores$Domain <- rownames(species_scores)
        
        # Create combined CCA plot (with environmental factors)
        all_plot <- ggplot() +
          # Add sample points
          geom_point(data = as.data.frame(all_scores$sites), 
                     aes(x = CCA1, y = CCA2), 
                     size = 3, color = "#93C3C6", alpha = 0.7) +
          
          # Add environmental factor arrows
          geom_segment(data = as.data.frame(all_scores$biplot), 
                       aes(x = 0, y = 0, xend = CCA1 * 2, yend = CCA2 * 2), 
                       arrow = arrow(length = unit(0.2, "cm")), 
                       color = "#E09D94") +
          
          # Add environmental factor labels
          geom_text(data = as.data.frame(all_scores$biplot), 
                    aes(x = CCA1 * 2.2, y = CCA2 * 2.2, label = rownames(all_scores$biplot)), 
                    color = "#E09D94", size = 4) +
          
          # Add microbial domain points (colored by domain)
          geom_point(data = species_scores, 
                     aes(x = CCA1, y = CCA2, color = Domain), 
                     size = 5, shape = 16) +
          
          # Add microbial domain labels
          geom_text(data = species_scores, 
                    aes(x = CCA1, y = CCA2, label = Domain), 
                    color = "black", size = 4, vjust = -1.2) +
          
          # Set plot titles and labels
          labs(title = paste0("All Microbes CCA Analysis (Adj R²: ", 
                              round(all_variance_explained, 2), "%)"),
               x = paste0("CCA1 (", round(all_cca$CCA$eig[1]/sum(all_cca$CCA$eig)*100, 2), "%)"),
               y = paste0("CCA2 (", round(all_cca$CCA$eig[2]/sum(all_cca$CCA$eig)*100, 2), "%)"),
               color = "Microbial Domain") +
          
          # Theme settings
          theme_minimal() +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
                axis.title = element_text(size = 12),
                legend.position = "bottom",
                legend.title = element_text(face = "bold"))
      }
      
      # Display combined CCA plot
      print(all_plot)
      
      # Export combined CCA results
      export_cca_results_csv(list(model = all_cca), soil_for_cca, "All Microbes", "all_microbes_cca_results")
    }
  } else {
    cat("Error: CCA model failed to fit properly, check data quality and variable selection.\n")
  }
}