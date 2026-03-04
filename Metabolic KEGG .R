
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

library('ComplexHeatmap')
library('circlize')
library('openxlsx')

# Read data
data <- read.xlsx("", sheet = 1, check.names = FALSE)

# Rename columns (replace APMS with NPSS)
colnames(data) <- gsub("APMS", "NPSS", colnames(data))

# Data processing
row.names(data) <- data$Name
data <- data[, -1]
data_matrix <- as.matrix(data)
scaled_data <- t(scale(t(data_matrix)))

# Define color gradient
color_gradient <- colorRamp2(c(-1.45, 0, 2.27), c("#0da9ce","white","#e74a32"))

# Create heatmap object
heatmap_obj <- Heatmap(
  scaled_data,
  col = color_gradient,
  name = "legend",
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 6, fontfamily = "serif"),
  column_names_gp = gpar(fontsize = 8, fontfamily ="serif"),
  heatmap_legend_param = list(
    title_gp = gpar(fontfamily = "serif", fontsize = 10),
    labels_gp = gpar(fontfamily = "serif", fontsize = 8)
  )
)

print(heatmap_obj)

# Set output directory
output_directory <- ""
if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)

# Save heatmap as PNG
png(file.path(output_directory, ""), 
    width = 13, height = 10, units = "cm", res = 1200)
draw(heatmap_obj)
dev.off()

# Save as PDF
pdf(file.path(output_directory, ""), 
    width = 13/2.54, height = 10/2.54)
draw(heatmap_obj)
dev.off()

# Save as TIFF
tiff(file.path(output_directory, ""), 
     width = 13, height = 10, units = "cm", res = 1200, compression = "lzw")
draw(heatmap_obj)
dev.off()