# Affiliation: Guozhongyuan Lab
# Correspondence: Guozhongyuan, guozhongyuan@sxmu.edu.cn

# Install required packages (uncomment first run)
# install.packages(c("tidyverse","ggrepel","openxlsx")) 
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ropls")

# Load libraries
library(openxlsx)
library(tidyverse)
library(ggrepel)
library(ropls)

# Read metabolomics data
metabo_data <- read.xlsx("", sheet = 1, check.names = FALSE)

# Extract numerical matrix
value_matrix <- metabo_data %>% 
  select(starts_with("all_")) %>% 
  as.matrix()
rownames(value_matrix) <- metabo_data$ID

# Create group information
group_info <- data.frame(
  sample = colnames(value_matrix),
  group = rep(c("APMS", "PSS"), each = 6)
)

# PLS-DA analysis
plsda_model <- opls(
  x = t(value_matrix),
  y = group_info$group,
  predI = 2,
  orthoI = 0,
  scaleC = "standard"
)

# Extract VIP scores
vip_scores <- getVipVn(plsda_model)

# Differential analysis
apms_means <- rowMeans(value_matrix[, 1:6])
pss_means <- rowMeans(value_matrix[, 7:12])

# Create volcano plot data frame
volcano_data <- data.frame(
  ID = rownames(value_matrix),
  Vip = vip_scores,
  FC = apms_means / pss_means,
  log2FC = log2(apms_means / pss_means),
  pvalue = apply(value_matrix, 1, function(x) {
    t.test(x[1:6], x[7:12])$p.value
  }), 
  stringsAsFactors = FALSE
) %>% 
  mutate(
    logp = -log10(pvalue),
    type = case_when(
      FC >= 2 & pvalue <= 0.05 & Vip > 1 ~ "up",
      FC <= 0.5 & pvalue <= 0.05 & Vip > 1 ~ "down",
      TRUE ~ "insig"
    )
  ) %>% 
  left_join(metabo_data %>% select(ID, Name), by = "ID")

# Save differential analysis results
write.xlsx(volcano_data, file = "")

# Extract top 5 VIP significant metabolites
top5_vip <- volcano_data %>%
  filter(type != "insig") %>%
  arrange(desc(Vip)) %>%
  slice(1:5)

# Calculate up/down regulated metabolite counts
type_counts <- volcano_data %>%
  filter(type %in% c("up", "down")) %>%
  count(type) %>%
  mutate(label = paste0(type, ": ", n))

# Generate final volcano plot
final_plot <- ggplot(volcano_data, aes(log2FC, logp)) +
  geom_point(
    aes(fill = type, size = Vip),
    shape = 21,
    alpha = 0.7,
    color = "transparent",
    stroke = 0.3
  ) +
  geom_point(
    data = top5_vip,
    aes(fill = type),
    shape = 21,
    color = "black",
    size = 4.5,
    stroke = 0.8,
    show.legend = FALSE
  ) +
  scale_size_continuous(
    range = c(2, 6),
    guide = guide_legend(override.aes = list(color = "black"))
  ) +
  scale_fill_manual(
    values = c("up" = "#E64B35", "down" = "#4DBBD5", "insig" = "grey70"),
    name = "Status"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_text_repel(
    data = top5_vip,
    aes(label = metabo_data$Name[match(ID, metabo_data$ID)]),
    max.overlaps = Inf,
    size = 3,
    family = "serif"
  ) +
  labs(
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("P-value")),
    title = "PLS-DA Volcano Plot"
  ) +
  annotate(
    geom = "text",
    x = Inf, y = Inf,
    label = paste("Up:", type_counts$n[type_counts$type == "up"]),
    color = "#E64B35",
    hjust = 1.5, vjust = 2.0,
    size = 5,
    family = "serif"
  ) +
  annotate(
    geom = "text", 
    x = -Inf, y = Inf,
    label = paste("Down:", type_counts$n[type_counts$type == "down"]),
    color = "#4DBBD5",
    hjust = -0.1, vjust = 2.0,
    size = 4.5,
    family = "serif"
  ) +
  theme_bw(base_family = "serif") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(size = 1.5, color = "black"),
    legend.justification = c(1, 0.5),
    legend.background = element_rect(
      fill = "transparent",
      color = "transparent",
      linewidth = 0.5
    ),
    legend.spacing.y = unit(0.2, "cm"),
    text = element_text(family = "serif" , size = 10),
    axis.title = element_text(face = "bold", size = 10),
    plot.title = element_text(hjust = 0.5, size = 15)
  ) 

print(final_plot)

# Save plot
ggsave(final_plot, file = "", width = 13, height = 10, units = 'cm', dpi = 1200)
ggsave(final_plot, file = "", width = 13, height = 10, units = 'cm', dpi = 1200)
ggsave(final_plot, file = "", width = 13, height = 10, units = 'cm', dpi = 1200)