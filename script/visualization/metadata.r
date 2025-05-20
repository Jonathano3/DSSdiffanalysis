library(readxl)
library(ggplot2)
library(dplyr)
library(viridis)

# Load Excel file
data <- read_excel("N:/donnee_poule/data_blood_novo/WP1_blood_novo_1149hens.xlsx") 

# Boxplot of raw read counts
distrib <- ggplot(data = data, aes(y = reads)) +
  geom_boxplot(fill = "#00a9c9", color = "black") +
  labs(title = "Boxplot des lectures brute", y = "nombre de lecures") +
  theme_minimal()

# Save boxplot
ggsave("N:\\figure\\boxplot_lecture_brute.png", plot = distrib)

# Density plot of read distribution per plate
figure3 <- ggplot(data = data, aes(x = reads, color = plateLibraryRRBS, fill = plateLibraryRRBS)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  scale_color_viridis_d(option = "H") + 
  scale_fill_viridis_d(option = "H") +
  labs(
    title = "Density Distribution of Reads by Plate",
    x = "Reads",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14)
  ) 

# Save density plot
ggsave("N:\\figure\\density_distribution.png", plot = figure3)

# Compute total and mean reads per barcode
stat_class2 <- data %>%
  group_by(barcodeSequence) %>% 
  summarise(
    total_reads = sum(reads),
    individuals = n()
  ) %>%
  mutate(mean_read = total_reads / individuals)

# Assign numeric factor to barcodes and sort data
data$barcode_num <- as.numeric(factor(data$barcodeSequence))
data <- data[order(data$barcode_num), ]

# Define axis breaks for barcode boxplot
breaks <- sort(unique(c(1, seq(5, max(data$barcode_num), by = 5))))

# Boxplot of reads by barcode
figure4_boxplot <- ggplot(data, aes(x = factor(barcode_num), y = reads)) +
  geom_boxplot(fill = "#00a9c9", color = "black", outlier.colour = "red", outlier.shape = 20) +
  theme_minimal() +
  labs(
    title = "Distribution des lectures en fonction des barcodes",
    x = "Barcode",
    y = "lecture"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_x_discrete(
    breaks = as.character(breaks),
    labels = breaks
  )

# Save barcode boxplot
ggsave("N:\\figure\\read_barcode.png", plot = figure4_boxplot)