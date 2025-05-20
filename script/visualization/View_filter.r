data <- read.delim("N:/donnee_poule/filtre/statistique_NA.txt")
stat <- read.delim("N:/donnee_poule/filtre/stat.txt")

library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

stat <- stat[stat$BED_PROCESSES > 180000,]
stat <- stat %>% mutate(
  filtre_min = BED - BED_PROCESSES
)

p1<- ggplot(stat, aes(x = BED_PROCESSES)) +
  geom_density(fill = "#2a9d8f", color = "grey", alpha = 0.5) +
  geom_vline(aes(xintercept = 180000, linetype = "Seuil choisi (180 000 CpG)"), color = "red", size = 1.5) +
  geom_vline(aes(xintercept = 150000, linetype = "Seuil envisagé (150 000/250 000 CpG)"), color = "#e76f51", size = 1) +
  geom_vline(aes(xintercept = 250000, linetype = "Seuil envisagé (150 000/250 000 CpG)"), color = "#e76f51", size = 1) +
  
  scale_linetype_manual(values = c("Seuil choisi (180 000 CpG)" = "solid", 
                                   "Seuil envisagé (150 000/250 000 CpG)" = "dashed")) +
  scale_y_continuous(labels = label_number()) + 
  labs(title = "Densité du nombre de CpG dans les échantillons bed_processed",
       x = "CpG",
       y = "Densité",
       linetype = "Seuil") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = c(0.82, 0.9),
  )

ggsave("N:\\figure\\seuil_bed_processed.png", plot = p1)

data <- data.frame(
  pipeline = 7826168,
  min = 2190123,
  low = 2189261,
  variant = 2087858,
  manquante = 592391
  
)

data_long <- data %>%
  select(pipeline, min, low,variant,manquante) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Pourcentage_NA")


min = data$pipeline - data$min
low = data$min - data$low
variant = data$low - data$variant
manquante = data$variant - data$manquante
all_data = data$manquante

values <- c(min, low, variant, manquante, all_data)
labels <- c("Filtre profondeur","Filtre qualité", "Filtre variant", "Filtre valeurs manquantes", "Données restantes")
colors <- c("#FF9999", "yellow", "#99FF99","#66B3FF", "#FFCC99")
percentages <- round(values / sum(values) * 100, 2)
labels_with_percent <- paste(labels, "\n", percentages, "%", sep = "")

png("N:/figure/filtre/valeur_NA_camembert.png", width = 800, height = 600)

par(mar = c(5, 4, 6, 2))
pie(values, labels = labels_with_percent, col = colors,
    main = "", 
    border = "white", cex = 1.2, radius = 1) 

title(main = "Répartition des données de CpG (en CpG unique)", line = 4, cex.main = 1.5)

dev.off()


