library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)   
library(GenomicFeatures)  
library(tidyr)
library(gridExtra)
library(GenomeFeatures)  
library(UpSetR)  
library(plotly)  
library(GenomicDistributions)  
library(txdbmaker)
# === Import genome annotation and prepare TxDb ===
gtf <- rtracklayer::import("/work/project/geronimo/fasta/1_ANOT_GRCg7b.GeneEnrichedAtlasFromENS107AndNCBI106.gtf", format="GFF")
txdb <- txdbmaker::makeTxDbFromGFF("/work/project/geronimo/fasta/1_ANOT_GRCg7b.GeneEnrichedAtlasFromENS107AndNCBI106.gtf", format='gtf', organism="Gallus gallus")

# === Load theoretical and final CpG positions ===
cpg_theo <- read.table("/work/project/geronimo/WP1/Jonathan/final/data/list_pos_theor.txt", header=FALSE, sep='\t')
colnames(cpg_theo) <- c('chr','start','end')
gr_cpg_theo <- GRanges(seqnames=cpg_theo$chr, ranges=IRanges(start=cpg_theo$start, end=cpg_theo$end))

cpg <- read.table("/work/project/geronimo/WP1/Jonathan/final/data/list_cpg_final_sb2.txt", header=FALSE, sep=',')
colnames(cpg) <- c('chr','start','end')
gr_cpg <- GRanges(seqnames=cpg$chr, ranges=IRanges(start=cpg$start, end=cpg$end))

# === Find common CpGs ===
common_CpG <- merge(cpg_theo, cpg, by=c("chr","start","end"))
common_CpG_gr <- GRanges(seqnames=common_CpG$chr, strand="*", ranges=common_CpG$start)

# === Calculate distance to TSS ===
tss <- promoters(gtf, upstream=0, downstream=1)
tss$gene_id <- names(tss)
strand(tss) <- "*"
tss <- unique(tss)
peaks <- GenomicRanges::GRangesList(theoric=gr_cpg_theo, final=gr_cpg)
distance_to_feature <- calcFeatureDist(peaks, tss)
for(i in names(peaks)) { peaks[[i]]$distance <- distance_to_feature[[i]] }

# === Build genome features ===
features2 <- GenomeFeatures::build_genome_features(
  txdb,
  features = c("promoter", "UTR5", "UTR3", "exons", "introns"),
  parameters = list(
    promoter_ranges = list(upstream=2000, downstream=0),
    downstream_range = 1000
  )
)

# === Overlap CpGs with genome features ===
features_overlaps <- suppressWarnings(
  GenomeFeatures::genome_features_overlaps(
    peaks,
    features2,
    figs_path = "/work/project/geronimo/WP1/Jonathan/final/plot/upset.png",
    ignore.strand = TRUE
  )
)

# === Plot distribution of CpGs across genome features ===
df <- data.frame(feature=rep(c('promoter','UTR5','UTR3','exons','introns','other'),2),
                 type=c(rep('Final',6),rep('Theoretical',6)),
                 percent=c(30.796,6.962,1.622,21.453,38.442,0.725,12.613,2.421,2.489,15.332,64.545,2.600))
df$type <- factor(df$type,levels=c('Theoretical','Final'))
df$feature <- factor(df$feature,levels=c('promoter','UTR5','UTR3','exons','introns','other'))
df <- df[order(df$type,rev(df$feature)),]

# === Annotate CpGs with CpG island context ===
cpg_file <- "/work/project/geronimo/fasta/ggal7_cpgIsland.bed.gz"
cpg_islands <- read.table(cpg_file, skip=1)
colnames(cpg_islands) <- c('chr','start','stop','nb_cpg')
conversion <- fread("/work/project/geronimo/fasta/sequence_report_ggal7.tsv")
cpg_islands$CHR <- NA
for(i in 1:nrow(cpg_islands)){
  if(grepl('NC_',cpg_islands$chr[i])){
    conv <- conversion$`Chromosome name`[conversion$`RefSeq seq accession`==cpg_islands$chr[i]]
    cpg_islands$CHR[i] <- conv
  }else{
    conv <- conversion$`GenBank seq accession`[conversion$`RefSeq seq accession`==cpg_islands$chr[i]]
    cpg_islands$CHR[i] <- conv
  }
}

gr_island <- GRanges(seqnames = cpg_islands$CHR, ranges = IRanges(start = cpg_islands$start, end = cpg_islands$stop))
north_shores_gr <- flank(gr_island, width = 2000, start = TRUE)
south_shores_gr <- flank(gr_island, width = 2000, start = FALSE)
north_shelves_gr <- flank(north_shores_gr, width = 2000, start = TRUE)
south_shelves_gr <- flank(south_shores_gr, width = 2000, start = FALSE)

annotation_cpg_theo <- rep("Open sea", length(gr_cpg_theo))
hits_islands <- findOverlaps(gr_cpg_theo, gr_island)
annotation_cpg_theo[queryHits(hits_islands)] <- "CpG island"
hits_north_shores <- findOverlaps(gr_cpg_theo, north_shores_gr)
annotation_cpg_theo[queryHits(hits_north_shores)] <- "North shore"
hits_south_shores <- findOverlaps(gr_cpg_theo, south_shores_gr)
annotation_cpg_theo[queryHits(hits_south_shores)] <- "South shore"
hits_north_shelves <- findOverlaps(gr_cpg_theo, north_shelves_gr)
annotation_cpg_theo[queryHits(hits_north_shelves)] <- "North shelf"
hits_south_shelves <- findOverlaps(gr_cpg_theo, south_shelves_gr)
annotation_cpg_theo[queryHits(hits_south_shelves)] <- "South shelf"
mcols(gr_cpg_theo)$annotation_CpG <- annotation_cpg_theo
CpG_isl_theor <- as.data.frame(gr_cpg_theo)
cpg_island_theo <- as.data.frame(table(CpG_isl_theor$annotation_CpG))
colnames(cpg_island_theo) <- c('region','count')
cpg_island_theo$type <- 'Theoretical'
cpg_island_theo$perc <- cpg_island_theo$count / sum(cpg_island_theo$count)

annotation_cpg <- rep("Open sea", length(gr_cpg))
hits_islands <- findOverlaps(gr_cpg, gr_island)
annotation_cpg[queryHits(hits_islands)] <- "CpG island"
hits_north_shores <- findOverlaps(gr_cpg, north_shores_gr)
annotation_cpg[queryHits(hits_north_shores)] <- "North shore"
hits_south_shores <- findOverlaps(gr_cpg, south_shores_gr)
annotation_cpg[queryHits(hits_south_shores)] <- "South shore"
hits_north_shelves <- findOverlaps(gr_cpg, north_shelves_gr)
annotation_cpg[queryHits(hits_north_shelves)] <- "North shelf"
hits_south_shelves <- findOverlaps(gr_cpg, south_shelves_gr)
annotation_cpg[queryHits(hits_south_shelves)] <- "South shelf"
mcols(gr_cpg)$annotation_CpG <- annotation_cpg
CpG_isl <- as.data.frame(gr_cpg)
cpg_island <- as.data.frame(table(CpG_isl$annotation_CpG))
colnames(cpg_island) <- c('region','count')
cpg_island$type <- 'Final'
cpg_island$perc <- cpg_island$count / sum(cpg_island$count)

cpg_island_theo$region <- factor(cpg_island_theo$region, levels=c('CpG island','North shore','South shore','North shelf','South shelf','Open sea'))
cpg_island <- cpg_island[order(cpg_island$region),]
cpg_island_plot <- rbind(cpg_island_theo, cpg_island)
cpg_island_plot$type <- factor(cpg_island_plot$type, levels=c('Theoretical','Final'))
cpg_island_plot <- cpg_island_plot[order(cpg_island_plot$type, rev(cpg_island_plot$region)),]

png('/work/project/geronimo/WP1/Jonathan/final/plot/perc_island_sb2.png', width=1000, height=1500, res=200)
ggplot(cpg_island_plot) +
  geom_bar(aes(x=type, y=perc, fill=region), position='stack', stat='identity') +
  theme_bw() +
  geom_text(aes(x=type, y=perc, label=round(perc, digit=2)), position="stack", vjust=+1.1, size=4)
dev.off()

# === Load methylation data and compute group stats ===
DF_sub2_meth <- fread("/work/project/geronimo/WP1/Jonathan/final/meth/DF_sub2_meth.txt", header=TRUE)
cage <- grep("1$", colnames(DF_sub2_meth))
no_cage <- grep("0$", colnames(DF_sub2_meth))
w70 <- grep("-70w-", colnames(DF_sub2_meth))
w90 <- grep("-90w-", colnames(DF_sub2_meth))

DF_sub2_meth$sd_meth_cage <- apply(DF_sub2_meth[, ..cage], 1, sd, na.rm = TRUE)
DF_sub2_meth$sd_meth_no_cage <- apply(DF_sub2_meth[, ..no_cage], 1, sd, na.rm = TRUE)
DF_sub2_meth$sd_meth_w70 <- apply(DF_sub2_meth[, ..w70], 1, sd, na.rm = TRUE)
DF_sub2_meth$sd_meth_w90 <- apply(DF_sub2_meth[, ..w90], 1, sd, na.rm = TRUE)

DF_sub2_meth$mean_70w_meth <- apply(DF_sub2_meth[,..w70], 1, mean, na.rm = TRUE)
DF_sub2_meth$mean_w90_meth <- apply(DF_sub2_meth[, ..w90], 1, mean, na.rm = TRUE)
DF_sub2_meth$mean_cage_meth <- apply(DF_sub2_meth[,..cage], 1, mean, na.rm = TRUE)
DF_sub2_meth$mean_no_cage_meth <- apply(DF_sub2_meth[,..no_cage], 1, mean, na.rm = TRUE)

# === Plot methylation variability: 70w vs 90w ===
plot_90w_70w <- DF_sub2_meth %>%
  mutate(delta_70w_90w = abs(mean_70w_meth - mean_w90_meth)) %>%
  select(delta_70w_90w, sd_meth, mean_meth) %>%
  drop_na()

p3 <- ggplot(plot_90w_70w, aes(x = mean_meth, y = sd_meth, color = delta_70w_90w)) +
  geom_point(alpha = 0.5) + 
  labs(x = "Moyenne", y = "Écart-type", title = "Écart-type sur moyenne méthylation (W70 vs W90)", color = "Delta mean 70w-90w") +
  theme_minimal() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient(low = "lightgrey", high = "red")

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/S.D_delta_70w_90w.png", plot = p3, width = 8, height = 6)

# === Plot methylation variability: cage vs no cage ===
plot_no_cage_cage <- DF_sub2_meth %>%
  mutate(delta_no_cage_cage = abs(mean_no_cage_meth - mean_cage_meth)) %>%
  select(delta_no_cage_cage, sd_meth, mean_meth) %>%
  drop_na()

p4 <- ggplot(plot_no_cage_cage, aes(x = mean_meth, y = sd_meth, color = delta_no_cage_cage)) +
  geom_point(alpha = 0.5) +
  labs(x = "Moyenne", y = "Écart-type", title = "Écart-type sur moyenne méthylation (no cage vs cage)", color = "Delta mean cage-no_cage") +
  theme_minimal() +
  theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5)) +
  scale_color_gradient(low = "lightgrey", high = "red")

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/S.D_delta_no_cage_cage.png", plot = p4, width = 8, height = 6)