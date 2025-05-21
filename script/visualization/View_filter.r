library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

# -------- Load data --------
data <- read.delim("N:/donnee_poule/filtre/statistique_NA.txt")
stat <- read.delim("N:/donnee_poule/filtre/stat.txt")

# -------- Filter samples with enough processed CpGs --------
stat <- stat %>% 
  filter(BED_PROCESSES > 180000) %>%
  mutate(filtre_min = BED - BED_PROCESSES)

# -------- Density plot of BED_PROCESSES with thresholds --------
p1 <- ggplot(stat, aes(x = BED_PROCESSES)) +
  geom_density(fill = "#2a9d8f", color = "grey", alpha = 0.5) +
  geom_vline(aes(xintercept = 180000, linetype = "Selected threshold (180k CpGs)"), 
             color = "red", size = 1.5) +
  geom_vline(aes(xintercept = 150000, linetype = "Alternative threshold (150k/250k CpGs)"), 
             color = "#e76f51", size = 1) +
  geom_vline(aes(xintercept = 250000, linetype = "Alternative threshold (150k/250k CpGs)"), 
             color = "#e76f51", size = 1) +
  scale_linetype_manual(values = c(
    "Selected threshold (180k CpGs)" = "solid", 
    "Alternative threshold (150k/250k CpGs)" = "dashed")) +
  scale_y_continuous(labels = label_number()) +
  labs(
    title = "Density of CpG counts in BED_PROCESSES samples",
    x = "CpGs",
    y = "Density",
    linetype = "Threshold"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.82, 0.9))

ggsave("N:\\figure\\seuil_bed_processed.png", plot = p1)

# -------- Donut/pie chart for NA filtering steps --------
data <- data.frame(
  pipeline = 7826168,
  min = 2190123,
  low = 2189261,
  variant = 2087858,
  manquante = 592391
)

# Compute filtered values step by step
val_min <- data$pipeline - data$min
val_low <- data$min - data$low
val_variant <- data$low - data$variant
val_na <- data$variant - data$manquante
val_final <- data$manquante

# Combine all values and labels
values <- c(val_min, val_low, val_variant, val_na, val_final)
labels <- c("Depth filter", "Quality filter", "Variant filter", "Missing values filter", "Remaining data")
colors <- c("#FF9999", "yellow", "#99FF99", "#66B3FF", "#FFCC99")

# Compute and format percentages
percentages <- round(values / sum(values) * 100, 2)
labels_with_percent <- paste(labels, "\n", percentages, "%")

# Save pie chart as image
png("N:/figure/filtre/valeur_NA_camembert.png", width = 800, height = 600)
par(mar = c(5, 4, 6, 2))
pie(values, labels = labels_with_percent, col = colors,
    border = "white", cex = 1.2, radius = 1)
title(main = "Breakdown of CpG filtering steps", line = 4, cex.main = 1.5)
dev.off()
