# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

# ===================== 1. Environment Configuration =====================
cat("Current working directory: ", getwd(), "\n")

# Load packages
suppressPackageStartupMessages({
  library(microeco)
  library(ape)
  library(magrittr)
  library(igraph)
  library(rgexf)
  library(RColorBrewer)
  library(WGCNA)
  library(scales)
  library(png)       
  library(grid)      
  library(grDevices) 
  library(dplyr)     
})
options(warn = -1) 

# Global parameters (for Plot window preview)
par(family = "Arial", ask = FALSE) # ask=FALSE: auto switch Plot window images

# ===================== 2. Data Preprocessing =====================
otu_table <- read.csv("", row.names = 1, sep = ",", fileEncoding = "UTF-8")
taxonomy_table <- read.csv("", row.names = 1, sep = ",", fileEncoding = "UTF-8")

taxonomy_table %<>% tidy_taxonomy
taxonomy_table[taxonomy_table == "" | taxonomy_table == "unclassified" | taxonomy_table == "未分类"] <- "Unclassified"
colnames(taxonomy_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

mt <- microtable$new(otu_table = otu_table, tax_table = taxonomy_table) 
mt$tax_table %<>% subset(Kingdom %in% c("k__Archaea", "k__Bacteria"))
mt$filter_pollution(taxa= c("mitochondria", "chloroplast"))
mt$tidy_dataset()

set.seed(9527)
rare_size <- min(mt$sample_sums())
mt$rarefy_samples(sample.size = rare_size)
cat("✅ Sample rarefaction completed: each sample standardized to", rare_size, "sequences\n")

mt$save_table(dirpath = "non_plastic_group_archaea_Basic_files", sep = ",")

# ===================== 3. Co-occurrence Network Analysis =====================
set.seed(9527) 
t1 <- trans_network$new(
  dataset = mt,
  cor_method = "spearman",
  use_WGCNA_pearson_spearman = TRUE,
  filter_thres = 0.0001
)

# Extract original correlation matrix
cor_matrix <- t1$res_cor_p[[1]]  
p_matrix <- t1$res_cor_p[[2]]    

# Calculate network
t1$cal_network(
  COR_p_thres = 0.05,
  COR_p_adjust = "fdr",  
  COR_cut = 0.6,         
  COR_optimization = FALSE
)

# Module detection
if(ecount(t1$res_network) > 0){ 
  t1$cal_module(method = "cluster_fast_greedy")
  module_num <- length(unique(V(t1$res_network)$module))
  cat("✅ Module detection completed, identified", module_num, "modules\n")
}else{
  V(t1$res_network)$module <- "Module1"
  module_num <- 1
  cat("⚠️ No valid edges in network, assigned default module\n")
}

# Node/edge table processing
t1$get_node_table(node_roles = FALSE)
t1$get_edge_table()
node_names <- V(t1$res_network)$name
t1$res_node_table$module <- V(t1$res_network)$module[match(t1$res_node_table$name, node_names)]
t1$res_node_table$degree <- degree(t1$res_network)[match(t1$res_node_table$name, node_names)]
t1$res_node_table[is.na(t1$res_node_table)] <- "Module1"

# Extract edge correlation coefficients
edge_list <- as_edgelist(t1$res_network)
cor_values <- numeric(nrow(edge_list))
for(i in 1:nrow(edge_list)){
  node1 <- edge_list[i,1]
  node2 <- edge_list[i,2]
  cor_values[i] <- cor_matrix[node1, node2]  
}
cat("✅ Extracted correlation coefficients for", length(cor_values), "edges, first 5 values:", head(cor_values), "\n")

# Save network files
t1$save_network(filepath = "non_plastic_group_archaea_network.gexf")
write.csv(t1$res_edge_table, "non_plastic_group_archaea_res_edge_table.csv", row.names = FALSE)
write.csv(t1$res_node_table, "non_plastic_group_archaea_res_node_table.csv", row.names = FALSE)
cat("✅ Network properties saved: Nodes =", vcount(t1$res_network), ", Edges =", ecount(t1$res_network), "\n")

# ===================== 4. Visualization (Plot preview + Multi-format save) =====================
# Define network plotting function
plot_network_archaea <- function() {
  par(bg = "white", mar = c(5, 5, 4, 10), mgp = c(3, 1, 0))  
  V(t1$res_network)$color <- color_mod[V(t1$res_network)$module]
  V(t1$res_network)$size <- scales::rescale(degree(t1$res_network), to = c(3, 10))
  E(t1$res_network)$color <- ifelse(cor_values > 0, col_pos, col_neg)
  
  plot(t1$res_network,
       layout = layout_with_fr,
       edge.alpha = 0.7, edge.width = 0.6,
       vertex.label = NA,
       vertex.frame.color = "#333333", vertex.frame.width = 0.6,
       main = "Co-occurrence network of archaea in non-plastic group (Genus level)",
       main.cex = 1.6, main.font = 2, main.col = "#333333",
       xlab = "Dimension 1 of network layout",
       xlab.cex = 1.2,
       ylab = "Dimension 2 of network layout",
       ylab.cex = 1.2,
       col.lab = "#333333",
       asp = 1
  )
  
  # Legend
  legend_items <- c(names(color_mod), "Positive correlation", "Negative correlation")
  legend_cols <- c(color_mod, col_pos, col_neg)
  legend("right",
         legend = legend_items,
         fill = legend_cols,
         title = "Modules / Correlation type",
         title.cex = 1.1, title.font = 2, title.col = "#333333",
         cex = 0.9, bty = "n",
         xpd = TRUE, inset = c(-0.15, 0),
         seg.len = 2, ncol = 1
  )
}

# Define node degree analysis plotting function
plot_node_degree_archaea <- function() {
  par(bg = "white", mar = c(5, 5, 4, 6), mgp = c(3, 1, 0))  
  node_degree <- as.numeric(t1$res_node_table$degree)
  node_names <- substr(t1$res_node_table$name, 1, 6)
  node_x <- 1:length(node_degree)
  
  plot(node_x, node_degree, 
       type = "p", pch = 16, col = "#6AA84F", cex = 1.2,
       main = "Node degree analysis of archaea network in non-plastic group",
       main.cex = 1.5, main.font = 2, main.col = "#333333",
       xlab = "Node ID", xlab.cex = 1.2,
       ylab = "Node degree (Number of connections)", ylab.cex = 1.2,
       col.lab = "#333333", xaxs = "i", yaxs = "i",
       xlim = c(0, max(node_x)+1), ylim = c(0, max(node_degree)+1),
       cex.axis = 0.9, col.axis = "#333333"
  )
  text(node_x, node_degree + 0.1, 
       labels = node_names, cex = 0.7, col = "#555555", 
       srt = 45, pos = 3
  )
  grid(col = "#E0E0E0", lty = 2, lwd = 0.8)
  legend("topright",
         legend = "Node degree: Number of edges connected to a node",
         bty = "n", cex = 0.9, text.col = "#333333",
         title = "Definition", title.cex = 1.0, title.font = 2,
         inset = c(0.02, 0.02)
  )
}

# 4.1 Color configuration
unique_modules <- unique(t1$res_node_table$module)
morandi_cols <- c("#A85BCB", "#FD7850", "#F3CE52", "#E89F58", "#C8373B", "#822E5C", "#C040AA")
color_mod <- morandi_cols[1:length(unique_modules)]  
names(color_mod) <- unique_modules
col_pos <- "#FD6F48"  # Positive correlation color
col_neg <- "#1A86C8"  # Negative correlation color

# 4.2 Co-occurrence network plot (Preview + Multi-format save)
cat("🔍 Displaying co-occurrence network in Plot window...\n")
# Step 1: Real-time preview in Plot window
plot_network_archaea()
cat("✅ Co-occurrence network displayed in Plot window!\n")

# Step 2: Save as PNG
png("non_plastic_group_archaea_network_plot.png", width = 14, height = 10, units = "in", res = 300, type = "cairo")
plot_network_archaea()
dev.off()
cat("✅ Co-occurrence network saved as PNG!\n")

# Step 3: Save as TIFF (Lossless LZW compression for journal submission)
tiff("non_plastic_group_archaea_network_plot.tiff", 
     width = 14, height = 10, units = "in", res = 300, 
     compression = "lzw", type = "cairo")
plot_network_archaea()
dev.off()
cat("✅ Co-occurrence network saved as TIFF (lossless compression)!\n")

# Step 4: Save as PDF (Vector format, infinitely scalable)
pdf("non_plastic_group_archaea_network_plot.pdf", width = 14, height = 10)
plot_network_archaea()
dev.off()
cat("✅ Co-occurrence network saved as PDF (vector format)!\n")

# 4.3 Node degree analysis plot (Preview + Multi-format save)
cat("🔍 Displaying node degree analysis in Plot window...\n")
# Step 1: Real-time preview in Plot window
plot_node_degree_archaea()
cat("✅ Node degree analysis displayed in Plot window!\n")

# Step 2: Save as PNG
png("non_plastic_group_archaea_node_degree_plot.png", width = 12, height = 7, units = "in", res = 300, type = "cairo")
plot_node_degree_archaea()
dev.off()
cat("✅ Node degree analysis saved as PNG!\n")

# Step 3: Save as TIFF (Lossless LZW compression)
tiff("non_plastic_group_archaea_node_degree_plot.tiff", 
     width = 12, height = 7, units = "in", res = 300, 
     compression = "lzw", type = "cairo")
plot_node_degree_archaea()
dev.off()
cat("✅ Node degree analysis saved as TIFF (lossless compression)!\n")

# Step 4: Save as PDF (Vector format)
pdf("non_plastic_group_archaea_node_degree_plot.pdf", width = 12, height = 7)
plot_node_degree_archaea()
dev.off()
cat("✅ Node degree analysis saved as PDF (vector format)!\n")

# ===================== 5. Print Plot window image (Core function) =====================
cat("\n🖨️  Preparing to print Plot window image...\n")
# Print function (convert to PDF first for printer compatibility)
print_plot <- function(){
  tryCatch({
    dev.print(device = "pdf", file = "temp_print.pdf") 
    # Uncomment below for direct printing to default printer
    # dev.print() 
    cat("✅ Plot window image print command sent!\n")
  }, error = function(e){
    cat("❌ Print failed:", e$message, "\n")
    cat("💡 Solution: Check printer connection, or save as PNG and print manually\n")
  })
}
# Execute print (runs print for current Plot window image)
print_plot()

# ===================== 6. Optional: Display saved PNG in Plot window =====================
show_saved_png <- function(png_path){
  img <- readPNG(png_path)
  dev.new(width = dim(img)[2]/300, height = dim(img)[1]/300, unit = "in")
  grid.raster(img)
  cat(paste0("✅ Displayed in Plot window:", png_path, "\n"))
}
# Examples:
# show_saved_png("non_plastic_group_archaea_network_plot.png")
# show_saved_png("non_plastic_group_archaea_node_degree_plot.png")

# ===================== 7. Publication-level statistical summary =====================
network_stats <- data.frame(
  Metric = c(
    "Number of nodes (Genus level)",
    "Number of edges",
    "Number of modules",
    "Mean node degree",
    "Number of positive edges (ρ>0.6)",
    "Number of negative edges (ρ<-0.6)",
    "Correlation threshold",
    "P-value threshold (FDR adjusted)"
  ),
  Value = c(
    vcount(t1$res_network),
    ecount(t1$res_network),
    module_num,
    round(mean(as.numeric(t1$res_node_table$degree), na.rm = TRUE), 2),
    sum(cor_values > 0),
    sum(cor_values < 0),
    "|Spearman's ρ| ≥ 0.6",
    "FDR < 0.05"
  )
)
write.csv(network_stats, "non_plastic_group_archaea_network_stats_iMeta.csv", row.names = FALSE)
