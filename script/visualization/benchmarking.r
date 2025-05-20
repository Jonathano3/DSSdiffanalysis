library(ggplot2)
library(dplyr)
library(stringr)
library(viridis)
library(purrr)

# Load and label data from multiple individuals
data_a <- read.delim("N:/donnee_poule/data_blood_novo/B258-70w-blood-RRBS.biscuit.bed.gz", header=FALSE); data_a$V6 <- "B258"
data_b <- read.delim("N:/donnee_poule/data_blood_novo/B453-blood-70w-RRBS.biscuit.bed.gz", header=FALSE); data_b$V6 <- "B453"
data_c <- read.delim("N:/donnee_poule/data_blood_novo/C6-90w-blood-RRBS.biscuit.bed.gz", header=FALSE);  data_c$V6 <- "C6"
data_d <- read.delim("N:/donnee_poule/data_blood_novo/C92-sang-90w-RRBS.biscuit.bed.gz", header=FALSE); data_d$V6 <- "C92"
data_e <- read.delim("N:/donnee_poule/data_blood_novo/C150-blood-90w-RRBS.biscuit.bed.gz", header=FALSE); data_e$V6 <- "C150"
data_f <- read.delim("N:/donnee_poule/data_blood_novo/C127-blood-90w-RRBS.biscuit.bed.gz", header=FALSE); data_f$V6 <- "C127"

data <- rbind(data_a, data_b, data_c, data_d, data_e, data_f)

# Group special chromosome names
data <- data %>% mutate(V1 = case_when(
  str_starts(V1, "MU") ~ "MU",
  str_starts(V1, "JAEN") ~ "JAEN",
  TRUE ~ V1
))

# Split data per chromosome
data_list <- split(data, data$V1)

# -------- Threshold with second derivative --------
top_profondeur2 <- list()
for (i in seq_along(data_list)) {
  dens <- density(data_list[[i]]$V5)
  derive <- diff(dens$y) / diff(dens$x)
  derive_seconde_derive <- diff(derive) / diff(dens$x[-1])

  pente_seuil <- 0.0000005
  seuil_seconde <- -0.0012

  # Detect where the second derivative drops below the threshold
  if (any(derive_seconde_derive < seuil_seconde)) {
    valeur_seuil <- which(abs(derive) <= pente_seuil)[1]
  } else {
    valeur_seuil <- 999999999
  }

  # Save threshold value
  top_profondeur2[[i]] <- list(
    seuil = dens$x[valeur_seuil],
    nom = data_list[[i]][1, 1]
  )

  # Save diagnostic plot
  p1 <- plot(dens$x[-c(1, length(dens$x))], derive_seconde_derive, type = "l", col = "green", lwd = 2,
             main = paste("Second derivative", data_list[[i]][1, 1]))
  ggsave(paste0("N:\\figure\\profondeur\\derive_seconde\\plot_", i, ".png"), plot = p1)
}

# -------- Threshold with 0.1% quantile --------
top_profondeur <- list()
for (i in seq_along(data_list)) {
  quantile_seuil <- quantile(data_list[[i]]$V5, 0.999, na.rm = TRUE)
  top_profondeur[[i]] <- list(
    seuil = quantile_seuil,
    nom = data_list[[i]][1, 1]
  )
}

# -------- Plot functions --------

# Density plot by chromosome
create_plot_density <- function(data, seuil) {
  ggplot(data, aes(x = V5, fill = factor(V6))) +
    geom_density(alpha = 0.5) +
    scale_fill_viridis_d(option = "turbo") +
    scale_x_continuous(trans = "log1p") +
    geom_vline(aes(xintercept = seuil$seuil), linetype = "dashed", color = "red") + # threshold line
    labs(title = "Depth density with threshold", x = "Depth (log)", y = "Density") +
    facet_grid(V1 ~ ., scales = "free") +
    theme_minimal()
}

# Scatter plot of depth per position
create_plot <- function(data, seuil) {
  ggplot(data, aes(x = V2, y = V5, color = factor(V6))) +
    geom_line(alpha = 0.3) +
    geom_point(alpha = 0.3) +
    scale_color_viridis_d(option = "turbo") +
    geom_hline(aes(yintercept = seuil$seuil), linetype = "dashed", color = "red") + # threshold line
    labs(title = "Depth across positions", x = "Position", y = "Depth") +
    facet_grid(V1 ~ ., scales = "free_x") +
    theme_minimal()
}

# -------- Apply and save plots --------
# With second derivative thresholds
plots_density_derive <- mapply(create_plot_density, data_list, top_profondeur2, SIMPLIFY = FALSE)
plots_derive <- mapply(create_plot, data_list, top_profondeur2, SIMPLIFY = FALSE)

for (i in seq_along(plots_derive)) {
  chrom_name <- names(data_list)[i]
  ggsave(paste0("N:\\figure\\profondeur\\derive\\plot_density_", chrom_name, ".png"), plot = plots_density_derive[[i]])
  ggsave(paste0("N:\\figure\\profondeur\\derive\\plot_", chrom_name, ".png"), plot = plots_derive[[i]])
}

# With quantile thresholds
plots_density_quantile <- mapply(create_plot_density, data_list, top_profondeur, SIMPLIFY = FALSE)
plots_quantile <- mapply(create_plot, data_list, top_profondeur, SIMPLIFY = FALSE)

for (i in seq_along(plots_quantile)) {
  chrom_name <- names(data_list)[i]
  ggsave(paste0("N:\\figure\\profondeur\\quantile\\0.1%\\plot_density_quantile_0.1%_", chrom_name, ".png"), plot = plots_density_quantile[[i]])
  ggsave(paste0("N:\\figure\\profondeur\\quantile\\0.1%\\plot_0.1%_", chrom_name, ".png"), plot = plots_quantile[[i]])
}

## second analysis

# Directory containing .bed.gz CpG files
input <- "N:/donnee_poule/data_blood_novo"

# Read all .bed.gz files in the directory
fichiers <- list.files(path = input, pattern = "\\.bed\\.gz$", full.names = TRUE)

# Load and merge all data files
data_list <- lapply(fichiers, function(i) {
  df <- read.delim(i, header = FALSE)
  colnames(df)[1] <- "chrom"
  return(df)
})
data <- do.call(rbind, data_list)

# Merge special chromosome names into common labels and filter them out
data <- data %>%
  mutate(chrom = case_when(
    str_starts(chrom, "MU") ~ "MU",
    str_starts(chrom, "JAEN") ~ "JAEN",
    TRUE ~ chrom
  )) %>%
  filter(!chrom %in% c("MU", "JAEN"))  # Exclude unneeded chromosomes

# Count number of CpGs per chromosome
cpg_counts <- data %>%
  group_by(chrom) %>%
  summarise(nb_CpG = n(), .groups = "drop")

# Chromosome sizes from NCBI
taille <- c(196449156, 149539284, 110642502, 90861225, 59506338, 36220557, 36382834, 29578256,
            23733309, 20453248, 19638187, 20119077, 17905061, 15331188, 12703657, 2706039,
            11092391, 11623896, 10455293, 14265659, 6970754, 4686657, 6253421, 6478339,
            3067737, 5349051, 5228753, 5437364, 726478, 755666, 2457334, 125424, 3839931,
            3469343, 554126, 358375, 157853, 667312, 177356, 9109940, 86044486, 16784)

chromosome <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, "W", "Z", "MT")

ncbi <- data.frame(chrom = as.character(chromosome), taille = taille)

# Merge CpG counts with chromosome sizes and compute CpG/length ratio
cpg_ratio <- left_join(cpg_counts, ncbi, by = "chrom") %>%
  mutate(ratio_CpG = nb_CpG / taille)

# Sort chromosomes in natural order
chr_order <- mixedsort(cpg_ratio$chrom)
cpg_ratio$chrom <- factor(cpg_ratio$chrom, levels = chr_order)
CPG <- cpg_ratio[order(cpg_ratio$chrom), ]

# Plot CpG density per chromosome
ratio <- ggplot(CPG, aes(x = chrom, y = ratio_CpG)) +
  geom_bar(stat = "identity", fill = "#00a9c9") +
  labs(title = "Ratio de CpG par chromosome", x = "Chromosome", y = "CpG / taille (bp)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("N:/figure/filtre/ratio_CpG_taille.png", plot = ratio, width = 8, height = 6)