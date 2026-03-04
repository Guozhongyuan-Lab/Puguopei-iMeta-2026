# Developed by: [Puguopei], First Author
# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

rm(list=ls())

library(ggplot2)
library(ggExtra)
library(vegan)
library(ggthemes)
library(readxl)

# Read ASV abundance table (Excel format)
data_raw <- read_excel("", sheet = 1)  

rownames(data_raw) <- data_raw[[1]]  
data <- data_raw[, -1]  
data <- as.matrix(data)  

# Calculate relative abundance
data <- data / apply(data, 2, sum)
data <- t(data)  

# Read group table (Excel format)
group_raw <- read_excel("", sheet = 1)  

colnames(group_raw)[colnames(group_raw) == "ID"] <- "Sample_ID"  

rownames(group_raw) <- group_raw$Sample_ID
group <- group_raw  

# Calculate Bray-Curtis distance
bray <- vegdist(data, method = 'bray')
bray_matrix <- as.matrix(bray)
write.table(bray_matrix, "", sep = "\t", quote = FALSE)

# Perform PCoA analysis
pcoa <- cmdscale(bray, k = 3, eig = TRUE)  

# Extract and organize PCoA coordinates
pcoa_data <- data.frame(pcoa$point)  
pcoa_data$Sample_ID <- rownames(pcoa_data)  
names(pcoa_data)[1:3] <- c("PCoA1", "PCoA2", "PCoA3")  

# Calculate explained variance for each dimension
eig_total <- sum(pcoa$eig)  
eig_percent <- round((pcoa$eig / eig_total) * 100, 1)  

# Merge PCoA coordinates with group information
pcoa_result <- merge(pcoa_data, group, by = "Sample_ID")

# PERMANOVA analysis
dune.div <- adonis2(data ~ group, data = group, permutations = 999, method = "bray")
dune_adonis <- paste0("Adonis R²: ", round(dune.div$R2[1], 2), "; P-value: ", round(dune.div$`Pr(>F)`[1], 3))

# Draw basic PCoA scatter plot
p <- ggplot(pcoa_result, aes(x = PCoA1, y = PCoA2, color = group)) +
  geom_point(aes(shape = group), size = 5) +
  labs(x = paste("PCoA 1 (", eig_percent[1], "%)", sep = ""),
       y = paste("PCoA 2 (", eig_percent[2], "%)", sep = ""),
       caption = dune_adonis) +
  scale_colour_manual(values = c("#2f7eb0", "#d00000")) +
  theme(legend.position = c(0.9, 0.19),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(color = "black", fill = "transparent"),
        axis.text = element_text(color = "black", size = 10)) +
  geom_hline(aes(yintercept = 0), colour = "#BEBEBE", linetype = "dashed") +
  geom_vline(aes(xintercept = 0), colour = "#BEBEBE", linetype = "dashed")

# Save basic PCoA plot
ggsave("", plot = p, width = 5, height = 5, dpi = 600)  
ggsave("", plot = p, width = 5, height = 5, dpi = 600)  
ggsave("", plot = p, width = 5, height = 5)  

# Add confidence ellipse
p_with_ellipse <- p + stat_ellipse(data = pcoa_result,
                                   geom = "polygon",
                                   level = 0.9,
                                   linetype = 2,
                                   linewidth = 0.5,
                                   aes(fill = group),
                                   alpha = 0.3,
                                   show.legend = TRUE) +
  scale_fill_manual(values = c("#2f7eb0", "#d00000"))
print(p_with_ellipse)

# Save plot with confidence ellipse
ggsave("", plot = p_with_ellipse, width = 5, height = 5, dpi = 600)  
ggsave("", plot = p_with_ellipse, width = 5, height = 5, dpi = 600)  
ggsave("", plot = p_with_ellipse, width = 5, height = 5)  

# Add marginal density plot and save
png(file = "", width = 5, height = 5, res = 600, units = "in")  
p_marginal <- ggMarginal(
  p_with_ellipse,
  type = "density",
  margins = "both",
  size = 3.5,
  groupColour = FALSE,
  groupFill = TRUE
)
print(p_marginal)
dev.off()  

tiff(file = "", width = 5, height = 5, res = 600, units = "in")  
p_marginal <- ggMarginal(
  p_with_ellipse,
  type = "density",
  margins = "both",
  size = 3.5,
  groupColour = FALSE,
  groupFill = TRUE
)
print(p_marginal)
dev.off()  

pdf(file = "", width = 5, height = 5)  
p_marginal <- ggMarginal(
  p_with_ellipse,
  type = "density",
  margins = "both",
  size = 3.5,
  groupColour = FALSE,
  groupFill = TRUE
)
print(p_marginal)
dev.off()