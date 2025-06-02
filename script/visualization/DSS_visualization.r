# 1. Load required packages 
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)
library(reshape2)
library(clusterProfiler)
library(gtools)

# 2. Load DSS differential results
multifactor_int_cage <- fread(".../multifactor_int_cage.txt")
multifactor_int_w <- fread(".../multifactor_int_w.txt")
multifactor_int_interaction <- fread(".../multifactor_int_interaction.txt")
dmlTest_sm_cage_70 <- fread(".../ALL_dmlTest_sm_cage_70.txt")
dmlTest_sm_cage_90 <- fread(".../ALL_dmlTest_sm_cage_90.txt")

# 3. Merge results and compute metrics 
table <- data.table(
  chr = multifactor_int_cage$chr,
  pos = multifactor_int_cage$pos,
  fdr_cage = multifactor_int_cage$fdrs,
  fdr_w = multifactor_int_w$fdrs,
  fdr_interaction = multifactor_int_interaction$fdrs,
  fdr_cage_70 = dmlTest_sm_cage_70$fdr,
  fdr_cage_90 = dmlTest_sm_cage_90$fdr,
  sens_cage = ifelse(multifactor_int_cage$fdrs < 0.05, ifelse(multifactor_int_cage$stat < 0, "hypo", "hyper"),"non_signif"),
  sens_w = ifelse(multifactor_int_w$fdrs < 0.05, ifelse(multifactor_int_w$stat < 0, "hypo", "hyper"),"non_signif"),
  sens_interaction = ifelse(multifactor_int_interaction$fdrs < 0.05, ifelse(multifactor_int_interaction$stat < 0, "hypo", "hyper"),"non_signif"),
  sens_cage_70 = ifelse(dmlTest_sm_cage_70$fdr < 0.05, ifelse(dmlTest_sm_cage_70$diff < 0, "hypo", "hyper"),"non_signif"),
  sens_cage_90 = ifelse(dmlTest_sm_cage_90$fdr < 0.05, ifelse(dmlTest_sm_cage_90$diff < 0, "hypo", "hyper"),"non_signif")
) %>%
  mutate(chr = case_when(
    str_starts(chr, "MU") ~ "MU",
    str_starts(chr, "JAEN") ~ "JAEN",
    TRUE ~ chr
  ))

#  4. Create chromosome-wide positions for plotting 
chr_order <- mixedsort(unique(table$chr))
don <- table %>%
  mutate(chr = factor(chr, levels = chr_order)) %>%
  group_by(chr) %>%
  summarise(chr_len = max(pos), .groups = "drop") %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(table, ., by = "chr") %>%
  mutate(BPcum = pos + tot)

axisdf <- don %>%
  group_by(chr) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop")

#  5. Manhattan plot function 
make_manhattan <- function(df, fdr_col, sens_col, output_path, title = "") {
  gg <- ggplot(df, aes(x = BPcum, y = -log10(.data[[fdr_col]]))) +
    geom_point(aes(color = chr), alpha = 0.8, size = 1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), length(chr_order))) +
    geom_point(data = subset(df, .data[[sens_col]] == "hyper"),
               aes(x = BPcum, y = -log10(.data[[fdr_col]])), color = "red", size = 1.3) +
    geom_point(data = subset(df, .data[[sens_col]] == "hypo"),
               aes(x = BPcum, y = -log10(.data[[fdr_col]])), color = "blue", size = 1.3) +
    scale_x_continuous(label = axisdf$chr, breaks = axisdf$center) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(legend.position = "none", panel.border = element_blank(),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size = 10)) +
    ggtitle(title)
  ggsave(output_path, plot = gg, width = 100, height = 16, units = "cm", dpi = 300)
}

# 6. Generate Manhattan plots
make_manhattan(don, "fdr_cage", "sens_cage", ".../man_interaction_cage_sign.png")
make_manhattan(don, "fdr_w", "sens_w", ".../man_interaction_w_sign.png")
make_manhattan(don, "fdr_interaction", "sens_interaction", ".../man_interaction_interaction_sign.png")
make_manhattan(don, "fdr_cage_70", "sens_cage_70", ".../man_interaction_cage_70_sign.png")
make_manhattan(don, "fdr_cage_90", "sens_cage_90", ".../man_interaction_cage_90_sign.png")


# 7. Annotation step


# Load DML data to annotate
working <- multifactor_int_cage # you can change this depending on your input that u want

# Convert to GRanges
gr_dml <- GRanges(
  seqnames = working$chr,
  ranges = IRanges(start = working$pos, end = working$pos),
  pval = working$fdrs
)

# Prepare gene identifiers
DF_id <- genome_ref %>%
  select(id = ENSEMBL, name = NAME) %>%
  distinct(name, .keep_all = TRUE)

# Filter gene entries from the GFF file based on gene names
gff_filtered_2 <- gff[mcols(gff)$type == "gene"]
gff_filtered_2_gene <- gff_filtered_2[mcols(gff_filtered_2)$gene_name %in% DF_id$name]

# Add Ensembl gene IDs to GFF entries
mcols(gff_filtered_2_gene)$gene_id <- sapply(
  mcols(gff_filtered_2_gene)$gene_name,
  function(x) DF_id$id[DF_id$name == x]
)

# Create promoter regions (2 kb upstream of the TSS)
gff_filtered_2_promoter <- gff_filtered_2_gene
plus_strand <- strand(gff_filtered_2_promoter) == "+"
minus_strand <- strand(gff_filtered_2_promoter) == "-"

end(gff_filtered_2_promoter)[plus_strand] <- start(gff_filtered_2_promoter)[plus_strand]
start(gff_filtered_2_promoter)[plus_strand] <- start(gff_filtered_2_promoter)[plus_strand] - 2000

start(gff_filtered_2_promoter)[minus_strand] <- end(gff_filtered_2_promoter)[minus_strand]
end(gff_filtered_2_promoter)[minus_strand] <- end(gff_filtered_2_promoter)[minus_strand] + 2000

# Find overlaps between DMLs and gene/promoter regions
hits_g <- findOverlaps(gr_dml, gff_filtered_2_gene)
hits_p <- findOverlaps(gr_dml, gff_filtered_2_promoter)

# Build annotated DML tables
annotated_dml_gene <- as_tibble(gr_dml[queryHits(hits_g)]) %>%
  mutate(gene_id = mcols(gff_filtered_2_gene)$gene_id[subjectHits(hits_g)]) %>%
  filter(pval < 0.05)

annotated_dml_promoter <- as_tibble(gr_dml[queryHits(hits_p)]) %>%
  mutate(gene_id = mcols(gff_filtered_2_promoter)$gene_id[subjectHits(hits_p)]) %>%
  filter(pval < 0.05)

#8. GO Enrichment analysis

# Prepare GO term tables
gene2go <- genome_ref %>% select(GO, ENSEMBL)
gene2name <- genome_ref %>% select(GO, FUNCTION)

# Enrichment in promoters
ego_p <- enricher(
  gene = annotated_dml_promoter$gene_id,
  universe = unique(genome_ref$ENSEMBL),
  TERM2GENE = gene2go,
  TERM2NAME = gene2name,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# Enrichment in gene bodies
ego_g <- enricher(
  gene = annotated_dml_gene$gene_id,
  universe = unique(genome_ref$ENSEMBL),
  TERM2GENE = gene2go,
  TERM2NAME = gene2name,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# 9. Visualisation

# GO term dot plots
dot_p <- dotplot(ego_p, showCategory = 10) + ggtitle("DML in promoter")
dot_g <- dotplot(ego_g, showCategory = 10) + ggtitle("DML in gene body")

# GO network plots
network_p <- cnetplot(ego_p, circular = TRUE, node_label = "all", cex_label_gene = 0.8)
network_g <- cnetplot(ego_g, circular = TRUE, node_label = "all", cex_label_gene = 0.8)

# Save plots
ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/final_dotplot_hyper_p.png", plot = dot_p, width = 10, height = 8)
ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/final_dotplot_hyper_g.png", plot = dot_g, width = 10, height = 8)
ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/final_network_hypo_p.png", plot = network_p, width = 8, height = 6)
ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/final_network_hypo_g.png", plot = network_g, width = 8, height = 6)
