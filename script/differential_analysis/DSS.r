# Load required libraries
library(GenomicRanges)
library(data.table)
library(reshape2)
library(AnnotationDbi)
library(rtracklayer) 
library(biomaRt)
library(dplyr)
library(DSS)
library(tidyverse)
library(readxl)
library(ggplot2)
library(clusterProfiler)
library(gtools)

# 1. Load methylation data and metadata
DF_meth <- fread("path/to/your/working/directory/data/all_meth_final_unique_CpG_sb2.txt", header = TRUE)
DF_depth <- fread("path/to/your/working/directory/data/all_depth_final_unique_CpG_sb2.txt", header = TRUE)
meta <- read_excel("path/to/your/working/directory/data/WP1_blood_novo_1149hens.xlsx")

# 2. Correct column names for mislabeled samples (you can skip this step)
colnames(DF_meth)[849] <- "C370-blood-RRBS-plate13-90w-1"
colnames(DF_meth)[856] <- "C377-blood-RRBS-plate13-90w-1"
colnames(DF_depth)[849] <- "C370-blood-RRBS-plate13-90w-1"
colnames(DF_depth)[856] <- "C377-blood-RRBS-plate13-90w-1"

# 3. Format methylation and depth data into long format and merge them
DF_depth_long <- DF_depth %>% pivot_longer(cols = -rs, names_to = "Sample", values_to = "Depth")
DF_meth_long <- DF_meth %>% pivot_longer(cols = -rs, names_to = "Sample", values_to = "Meth")
merged_df <- inner_join(DF_depth_long, DF_meth_long, by = c("rs", "Sample"))

chrom_merged <- colsplit(merged_df$rs, ':', c('chr', 'pos'))
start_merged <- colsplit(chrom_merged$pos, '_', c('start', 'end'))
merged_df <- merged_df %>%
  mutate(
    chr = chrom_merged$chr,
    pos = start_merged$start,
    X = round(Depth * Meth),
    X = ifelse(is.na(X), 0, X),
    N = Depth,
    N = ifelse(is.na(N), 0, N)
  ) %>% select(chr, pos, Sample, N, X) %>% arrange(chr, as.numeric(pos))

DF <- merged_df %>% split(.$Sample) %>%
  lapply(function(x) {
    x %>% select(-Sample) %>% as.data.table() %>%
      .[, pos := as.numeric(pos)] %>% .[, chr := as.character(chr)]
  })

col_names <- names(DF)

# ─────────────────────────────────────────────
# SIMPLE MODEL SECTION
# ─────────────────────────────────────────────

# 4. Prepare data for simple model comparison between 70w and 90w groups
DF_70 <- DF[grep("70w-", names(DF))]
DF_90 <- DF[grep("90w-", names(DF))]
Bsobj_70 = makeBSseqData(DF_70, c(names(DF_70)))
Bsobj_90 = makeBSseqData(DF_90, c(names(DF_90)))

# 5. Define sample groups for simple model
DF_names <- names(DF)
DF_name_70 <- grep("70w-", DF_names, value = TRUE)
DF_name_90 <- grep("90w-", DF_names, value = TRUE)
cage_70 <- grep("1$", DF_name_70, value = TRUE)
no_cage_70 <- grep("0$", DF_name_70, value = TRUE)
cage_90 <- grep("1$", DF_name_90, value = TRUE)
no_cage_90 <- grep("0$", DF_name_90, value = TRUE)

# 6. Run DML test for 70w and 90w groups
dmlTest_sm_cage_70 <- DMLtest(Bsobj_70, group1 = c(no_cage_70), group2 = c(cage_70), smoothing = TRUE)
dmlTest_sm_cage_90 <- DMLtest(Bsobj_90, group1 = c(no_cage_90), group2 = c(cage_90), smoothing = TRUE)

# 7. Save DML results from the simple model
fwrite(dmlTest_sm_cage_70, 'path/to/your/working/directory/DSS/ALL_dmlTest_sm_cage_70.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(dmlTest_sm_cage_90, 'path/to/your/working/directory/DSS/ALL_dmlTest_sm_cage_90.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# ─────────────────────────────────────────────
# MULTIFACTOR MODEL SECTION
# ─────────────────────────────────────────────

# 8. Identify matched samples for multifactor paired analysis
pairs <- sub("^(.*?-[^-]*-[^-]*-[^-]*)-.*", "\\1", col_names)
animal_ids <- meta$animalId[match(pairs, meta$sampleIdRRBS)]
duplicated_ids <- animal_ids[duplicated(animal_ids)]
duplicated_ids <- duplicated_ids[!is.na(duplicated_ids) & duplicated_ids != "NA"]
pair_indices <- which(animal_ids %in% duplicated_ids)
DF_pair <- DF[pair_indices]


# 9. Build the design matrix with cage, age, and individual pair
design_int <- data.table()
for (i in 1:length(DF_pair)) {
  nom <- names(DF_pair[i])
  cage <- substr(nom, nchar(nom), nchar(nom))
  w <- substr(nom, nchar(nom)-4, nchar(nom)-2)

  design_int <- rbindlist(list(design_int, data.table(
    cage = factor(cage),
    w = factor(w)
  )), use.names = TRUE)
}

# 10. Fit multifactor DSS model with interaction terms
Bsobj_int = makeBSseqData(DF_pair, c(names(DF_pair)))
DMLfit_int = DMLfit.multiFactor(Bsobj_int, design = design_int, formula = ~cage + w + cage:w, smoothing = TRUE)

# 11. Perform DML tests for each model term
multifactor_int_cage <- DMLtest.multiFactor(DMLfit_int, term = "cage")
multifactor_int_w <- DMLtest.multiFactor(DMLfit_int, term = "w")
multifactor_int_interaction <- DMLtest.multiFactor(DMLfit_int, term = "cage:w")

# 12. Save DML test results from multifactor model
fwrite(multifactor_int_cage, "path/to/your/working/directory/DSS/DML_multifactor_int_cage.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(multifactor_int_w, "path/to/your/working/directory/DSS/DML_multifactor_int_w.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(multifactor_int_interaction, "path/to/your/working/directory/DSS/DML_multifactor_int_interaction.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# 13. Format DML results for DMR detection
multifactor_int_cage <- multifactor_int_cage[, c(1:3, 5, 4)]
colnames(multifactor_int_cage) <- c("chr", "pos", "stat", "pval", "fdrs")
multifactor_int_w <- multifactor_int_w[, c(1:3, 5, 4)]
colnames(multifactor_int_w) <- c("chr", "pos", "stat", "pval", "fdrs")
multifactor_int_interaction <- multifactor_int_interaction[, c(1:3, 5, 4)]
colnames(multifactor_int_interaction) <- c("chr", "pos", "stat", "pval", "fdrs")

# 14. Call DMRs for each model term
DMR_multifactor_int_cage <- callDMR(multifactor_int_cage, p.threshold = 0.05)
DMR_multifactor_int_w <- callDMR(multifactor_int_w, p.threshold = 0.05)
DMR_multifactor_int_interaction <- callDMR(multifactor_int_interaction, p.threshold = 0.05)

# 15. Save DMR results
fwrite(DMR_multifactor_int_cage, "path/to/your/working/directory/DSS/DMR_multifactor_int_cage.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(DMR_multifactor_int_w, "path/to/your/working/directory/DSS/DMR_multifactor_int_w.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
fwrite(DMR_multifactor_int_interaction, "path/to/your/working/directory/DSS/DMR_multifactor_int_interaction.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")