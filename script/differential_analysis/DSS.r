library(data.table)
library(reshape2)
library(dplyr)
library(DSS)
library(ggplot2)
library(ggExtra)
library(GenomicRanges)
library(rtracklayer)
library(ChIPpeakAnno)
library(tidyverse)
library(ggVennDiagram)

#load data
DF_meth <- fread("/work/project/geronimo/WP1/Jonathan/final/meth/DF_sub2_meth.txt", header = TRUE)
DF_depth <- fread("/work/project/geronimo/WP1/Jonathan/final/data/all_depth_final_unique_CpG_sb2.txt", header = TRUE)

#NA col 850 et 857 de DF_meth (provisoire)
colnames(DF_meth)[320] <- "C370-blood-RRBS-plate13-90w-1"
colnames(DF_meth)[327]<- "C377-blood-RRBS-plate13-90w-1"
colnames(DF_depth)[320] <- "C370-blood-RRBS-plate13-90w-1"
colnames(DF_depth)[327]<- "C377-blood-RRBS-plate13-90w-1"

#match start + 1
chrom_meth <- colsplit(DF_meth$rs, ':',c('chr','pos'))
final_meth <- colsplit(chrom_meth$pos, '_',c('start','end'))
final_meth$start<- final_meth$start + 1
DF_meth$rs <- paste0(chrom_meth$chr,':',final_meth$start,'_',final_meth$end)
DF_meth$pos <- final_meth$start
DF_meth$chr <- chrom_meth$chr

## filtre NA
#na_rate <- colSums(is.na(DF_meth[,2:(ncol(DF_meth)-3)]))
#na_rate_filtered <- na_rate[na_rate<0.7*nrow(DF_meth)] 
#DF_meth <- DF_meth %>% select(c("rs", names(na_rate_filtered), "chr", "pos"))

#data formatting
#DF_depth <- DF_depth %>% select(c("rs", names(na_rate_filtered)))
#DF_depth$pos <- colsplit(DF_depth$rs, ':',c('chr','pos'))$chr
#DF_depth <- DF_depth[DF_depth$pos == 1, ]
#DF_depth <- DF_depth[, 1:(ncol(DF_depth) - 1)]
#DF_meth <- DF_meth[DF_meth$chr == 1, ]
#DF_meth <- DF_meth[, 1:(ncol(DF_meth) - 2)]

DF_meth = DF_meth[,1:(ncol(DF_meth)-5)]
DF <- list()
DF_depth_long <- DF_depth %>% pivot_longer(cols = -rs, names_to = "Sample", values_to = "Depth")
DF_meth_long <- DF_meth %>% pivot_longer(cols = -rs, names_to = "Sample", values_to = "Meth")
merged_df <- inner_join(DF_depth_long, DF_meth_long, by = c("rs", "Sample"))
chrom_merged <-colsplit(merged_df$rs, ':',c('chr','pos'))
start_merged <- colsplit(chrom_merged$pos, '_',c('start','end'))

merged_df <- merged_df %>%
  mutate(
    chr = chrom_merged$chr,
    pos = start_merged$start,
    X = round(Depth * Meth),
    X = ifelse(is.na(X), 0, X),
    N = Depth,
    N = ifelse(is.na(N), 0, N)
  ) %>% select(chr, pos, Sample, N, X) %>% arrange(chr, as.numeric(pos))

DF <- merged_df %>%split(.$Sample) %>%lapply(function(x) {x %>% select(-Sample) %>% as.data.table()%>%.[, pos := as.numeric(pos)] %>% .[, chr := as.character(chr)]})


#################################

#converting into bsseq object
BSobj = makeBSseqData( DF, c(names(DF)))

#creation of the design matrix
design <- data.table()
for (i in 1:length(DF)){
     nom <- names(DF[i])
     cage <- substr(nom, nchar(nom), nchar(nom))
     w <- substr(nom, nchar(nom)-4, nchar(nom)-2)
     pair <- sub("^[A-Za-z](.*?)-.*", "\\1", nom)  
     design <- rbind(design, data.table(cage = cage, w = w, pair = pair))
}

#### MODELE LINEAIRE

#information recovery
DF_names <- names(DF)
cage<-grep("1$", DF_names, value = TRUE)
no_cage<-grep("0$", DF_names, value = TRUE)
w70<-grep("70w-", DF_names, value = TRUE)
w90<-grep("90w-", DF_names, value = TRUE)

######DMC
dmlTest_sm_cage <- DMLtest(BSobj, group1=c(cage), group2=c(no_cage) ,smoothing=TRUE)

#table creation 
res_DMC <- dmlTest_sm_cage
res_DMC$padj <- p.adjust(res_DMC$pval, method = "BH")
res_DMC$start <- res_DMC$pos 
res_DMC$end <- res_DMC$pos + 1
res_DMC$pos <- res_DMC$start
res_DMC$pool1 <- ifelse(res_DMC$diff > 0, "UP", "DOWN")
res_DMC <- res_DMC[, c("chr", "start", "end", "mu1", "mu2", "diff", "diff.se", "pval", "padj", "pool1","pos")]
res_DMC$strand <- "+"

annotatedDMC <- GRanges(seqnames = res_DMC$chr, ranges = IRanges(res_DMC$end, res_DMC$end), strand = res_DMC$strand, name = res_DMC$pool1, score = res_DMC$diff)

#save
#fwrite(res_DMC,'/work/project/geronimo/WP1/Jonathan/final/DSS/res_DMC_cage.txt',col.names=T,row.names=F,quote=F, sep = "\t")

######DMR
res_DMR <- callDMR(dmlTest_sm_cage, p.threshold=0.01)
colnames(res_DMR) <- c("chr", "start", "end", "length", "nCG","meanMethy1", "meanMethy2", "diff.Methy", "areaStat", "pool1")
res_DMR <- as.data.frame(res_DMR)

#table creation 
res_DMR$pool1 <- ifelse(res_DMR$diff.Methy > 0, "UP", "DOWN")
res_DMR_grange <- GRanges(seqnames = res_DMR$chr, ranges = IRanges(res_DMR$start, res_DMR$end),name = res_DMR$pool1, score = res_DMR$diff)
gff <- import.gff3("/work/project/geronimo/WP1/Jonathan/useful/annotation_promoteur_fixed.gff")
res_DMR_grange_annot <- GRanges(seqnames = res_DMR$chr, ranges = IRanges(res_DMR$start, res_DMR$end))
gene_grange <- GRanges(seqnames = seqnames(gff)[gff@elementMetadata@listData$type == "gene"],IRanges(start(gff)[gff@elementMetadata@listData$type == "gene"],end = end(gff)[gff@elementMetadata@listData$type == "gene"],names = gff@elementMetadata@listData$Name[gff@elementMetadata@listData$type == "gene"]))
annotatedDMR <- suppressWarnings(annotatePeakInBatch(res_DMR_grange_annot, AnnotationData=gene_grange, output = "both"))




#save 
#fwrite(res_DMR,'/work/project/geronimo/WP1/Jonathan/final/DSS/res_DMR_cage.txt',col.names=T,row.names=F,quote=F, sep = "\t")
#fwrite(annotatedDMR,'/work/project/geronimo/WP1/Jonathan/final/DSS/annotatedDMR_cage.txt',col.names=T,row.names=F,quote=F, sep = "\t")


####################GRAPH####################

#### volcano DML


res_DMC$FC <- log2((exp(res_DMC$mu1) ) / (exp(res_DMC$mu2) ) )

###graph distibution FC

hist_plot <- ggplot(res_DMC, aes(x = FC)) +
  geom_histogram(binwidth = 0.01, color = "black", alpha = 0.7, position = "identity") +
  labs(title = "Distribution des log2FC (Cage vs No Cage)", 
       x = "log2FC", 
       y = "Fréquence") +
  theme_minimal()


ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/log2FC_distribution_histogram.png", plot = hist_plot, width = 8, height = 6)

###graph distribution p-value

plot_pval <- ggplot(res_DMC, aes(x = padj)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution des p-values ajusted", x = "p-value ajusted", y = "Fréquence") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/pval_distribution_histogram.png", plot = plot_pval, width = 8, height = 6)

###graph NA rate
na_rate_long <- na_rate %>% pivot_longer(everything(), names_to = "Sample", values_to = "NA_rate")
hist <- ggplot(na_rate_long, aes(x = NA_rate)) +
  geom_histogram(binwidth = 5, fill = "lightblue", color = "black", alpha = 0.7) +
  labs(title = "Histogramme des pourcentages de NA", x = "Pourcentage de NA", y = "Fréquence") +
  theme_minimal()

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/na_rate_histogram.png", plot = hist, width = 8, height = 6)

###graph volcano

hist_vol <- ggplot(res_DMC, aes(x = FC, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), size = 1, alpha = 0.7) +
  scale_color_manual(values = c("black", "red")) + 
  labs(title = "Volcano plot", x = "FC", y = "-log10(padj)") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/volcano_cage_vs_nocage.png", plot = hist_vol, width = 8, height = 6)

###graph FC/diff
hist_dif <- ggplot(res_DMC2, aes(x = FC, y = diff)) +
  geom_point(alpha = 0.5, color = "blue") +
  labs(title = "FC vs diff", x = "FC", y = "diff") +
  theme_minimal()

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/FC_sur_diff.png", plot = hist_dif, width = 8, height = 6)

## qqplot

png("/work/project/geronimo/WP1/Jonathan/final/plot/qqplot_pvalues.png", 
    width = 1600, height = 1200, res = 300)
qqplot(runif(length(res_DMC2$pval)), res_DMC2$pval,
       main = "QQ-plot des p-values", 
       xlab = "Quantiles théoriques (Uniforme)", 
       ylab = "Quantiles observés (p-values)", 
       cex.main = 1, cex.lab = 0.9, cex.axis = 0.8)
abline(0, 1, col = "red", lwd = 2, lty = 2)

dev.off()



################################################################
#volcano plot

#data formatting
DF_volcano <- merge(res_DMC[ ,c("chr","pos","padj")], DT_count[, c("chr", "pos", "log2FC_cage")], by = c("chr", "pos"), all.x = TRUE)

#graph cage
plot <- ggplot(DF_volcano, aes(x = log2FC_cage, y = -log10(padj))) +
   geom_point(aes(color = padj < 0.05), size = 1, alpha = 0.7) + 
   scale_color_manual(values = c("black", "red")) +
   theme_minimal() +
   labs(title = "Volcano Plot des DML pour l'effet de la cage", x = "log2FC", y = "-log10(p-value)") +
   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")

volcano_marginal <- ggMarginal(plot, type = "density", fill = "gray", alpha = 0.5)

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/volcano_dml_cage.png", plot = volcano_marginal, width = 8, height = 6)

#graph w

DF_volcano <- merge(res_DMC[ ,c("chr","pos","padj","diff.se")], DT_count[, c("chr", "pos", "log2FC_w")], by = c("chr", "pos"), all.x = TRUE)

plot <- ggplot(DF_volcano, aes(x = log2FC_w, y = -log10(padj))) +
   geom_point(aes(color = diff.se > 0.001), size = 1, alpha = 0.7) + 
   scale_color_manual(values = c("black", "red")) + 
   theme_minimal() +
   labs(title = "Volcano Plot des DML pour l'effet de la duré", x = "log2FC", y = "-log10(p-value)") +
   geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray")

volcano_marginal <- ggMarginal(plot, type = "density", fill = "gray", alpha = 0.5)

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/volcano_dml_w.png", plot = volcano_marginal, width = 8, height = 6)


################################################################

###plot test
BSobj_smoothed <- BSmooth(BSobj)
png("/work/project/geronimo/WP1/Jonathan/final/plot/region_plot.png", width = 2800, height = 800, res = 300)
plotManyRegions(BSobj_smoothed, res_DMR[1:100,], extend = 5000)
dev.off()

##graph distance between 2 cytosines
diff <- unlist(sapply(unique(as.character(seqnames(BSobj))), function(x) diff(as.numeric(start(BSobj))[as.character(seqnames(BSobj)) == x[1]])))

p<-ggplot(data.frame(diff=diff), aes(x = log2(diff))) + 
    geom_density(adjust = 2) + ylab("Density") +
    theme_bw() + xlab(label = "Distance between 2 cytosines (log2)") + ggtitle("Distribution of the distance between 2 cytosines") 

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/distance_between_cytosine.png", plot = p, width = 8, height = 6)

###graph distance between DMR
res_DMC <- res_DMC[order(res_DMC$chr, res_DMC$start), ]
diffDMC <- unlist(sapply(as.character(unique(res_DMC$chr)), function(x) diff(res_DMC$start[res_DMC == x[1]])))
p2<-ggplot(data.frame(diffDMC=diffDMC), aes(x = log2(diffDMC))) + 
      geom_density() + ylab("Density") +
      theme_bw() + xlab(label = "Distance between 2 DMC (log2)") + ggtitle("Distribution of the distance between 2 DMC")

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/distance_between_DMC.png", plot = p2, width = 8, height = 6)

##numeber of DMC by type
dmc_type <- NULL
prop_by_type <- NULL
gff <- import.gff3("/work/project/geronimo/WP1/Jonathan/useful/annotation_promoteur_fixed.gff")
annot <- data.frame(chr = seqnames(gff), start = start(gff), end = end(gff), strand = strand(gff),type = gff@elementMetadata@listData$type)
type <- unlist(strsplit('gene,exon,five_prime_utr,three_prime_utr,promoter,repeat,TSS', ','))

annot <- annot[annot$type %in% type, ]
dmc_grange <- GRanges(seqnames = res_DMC$chr,IRanges(res_DMC$end, res_DMC$end),strand = res_DMC$strand)
dmc_type <- list()
for(i in type){
        annot_grange <- GRanges(seqnames = annot$chr[annot$type == i],
                                IRanges(annot$start[annot$type == i],
                                        annot$end[annot$type == i]),
                                strand = annot$strand[annot$type == i])
        subset_type  <- suppressWarnings(subsetByOverlaps(dmc_grange, annot_grange))
        dmc_type[[i]] <- paste0(seqnames(subset_type), ".", start(subset_type))
      }
dmc_type[["other"]] <- paste0(res_DMC$chr, ".", res_DMC$end)[!(paste0(res_DMC$chr, ".", res_DMC$end) %in% unlist(dmc_type))]
dmc_type <- dmc_type[dmc_type != "."]
number_type <- data.frame(type = as.factor(names(dmc_type)), frequency=sapply(dmc_type, length))
number_type$type <-factor(number_type$type, levels=number_type$type[order(number_type$frequency, decreasing = TRUE)])
p <- ggplot(number_type, aes(type, frequency)) + geom_bar(stat = "identity") +
        xlab("Type of region") + ylab("Frequency") + ggtitle("Number of DMC by type of region") +
        theme_bw()
ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/number_DMC_by_type.png", plot = p, width = 8, height = 6)

##venn diagram

dmc_type_filtered <- dmc_type[!names(dmc_type) %in% c( "gene","other")]

venn_plot <- ggVennDiagram(dmc_type_filtered, 
                   label_alpha = 0, 
                   stroke_size = 0.5,  
                   label_size = 3) +
  scale_fill_gradient(low = "#2a7fb2", high = "#f97373") + 
  ggtitle("Diagramme de Venn") +
  theme_minimal()

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/venn_diagram.png", plot = venn_plot, width = 8, height = 6, dpi = 300)

## distance DMC to TSS
annot <- data.frame(chr = seqnames(gff), start = start(gff), end = end(gff), strand = strand(gff),type = gff@elementMetadata@listData$type)
type <- unlist(strsplit('TSS', ','))
annot <- annot[annot$type %in% type, ]
annot <- annot[annot$chr == 1, ]
tssGrange <- GRanges(seqnames = annot$chr, IRanges(annot$start - 5000, annot$end + 5000))
res_DMC_sig <- res_DMC[res_DMC$padj < 0.05, ]
annotatedDMC <- GRanges(seqnames = res_DMC_sig$chr, ranges = IRanges(res_DMC_sig$end, res_DMC_sig$end))
overlaps_TSS <- suppressWarnings(subsetByOverlaps(tssGrange, annotatedDMC))
overlaps_DMC <- suppressWarnings(subsetByOverlaps(annotatedDMC, tssGrange))

adjusted_tss_start <- start(overlaps_TSS) + 5000
dmc_distance <- data.frame(start_dmc = numeric(length(overlaps_DMC)),distance_to_tss = numeric(length(overlaps_DMC)))
for (line in 1:length(overlaps_DMC)) {
  dmc_start <- start(overlaps_DMC)[line]
  distance <- dmc_start - adjusted_tss_start
  min_distance <- distance[which.min(abs(distance))]
  dmc_distance$start_dmc[line] <- dmc_start
  dmc_distance$distance_to_tss[line] <- min_distance
}

distance <- ggplot(dmc_distance, aes(x = distance_to_tss)) +
  geom_histogram(binwidth = 500, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution des distances entre DMC et TSS",
       x = "Distance au TSS ",
       y = "Nombre de DMC") +
  theme_minimal()

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/distance_DMC_TSS.png", plot = distance, width = 8, height = 6)


##############################################

##############################################
#save

fwrite(dmr_cage,'/work/project/geronimo/WP1/Jonathan/final/DSS/dmr_cage.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(dmr_w,'/work/project/geronimo/WP1/Jonathan/final/DSS/dmr_w.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(dmr_WxCAGE,'/work/project/geronimo/WP1/Jonathan/final/DSS/dmr_WxCAGE.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(DMLtest_cage,'/work/project/geronimo/WP1/Jonathan/final/DSS/DMLtest_cage.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(DMLtest_w,'/work/project/geronimo/WP1/Jonathan/final/DSS/DMLtest_w.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(DMLtest_WxCAGE,'/work/project/geronimo/WP1/Jonathan/final/DSS/DMLtest_WxCAGE.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(design,'/work/project/geronimo/WP1/Jonathan/final/data/design.txt',col.names=T,row.names=F,quote=F, sep = "\t")



###MODELE INTERACTION
# testing if effects have interaction between conditions
DMLfit_interaction = DMLfit.multiFactor(BSobj, design=design, formula=~cage+w+cage:w)


# testing effect of cage
DMLtest_cage = DMLtest.multiFactor(DMLfit_interaction, coef=2)
# testing effect of w
DMLtest_w = DMLtest.multiFactor(DMLfit_interaction, coef=3)
# testing interaction
DMLtest_WxCAGE = DMLtest.multiFactor(DMLfit_interaction, coef=4)
# DMR
dmr_cage = callDMR(DMLtest_cage, p.threshold=0.05)
dmr_w = callDMR(DMLtest_w, p.threshold=0.05)
dmr_WxCAGE = callDMR(DMLtest_WxCAGE, p.threshold=0.05)


########
library(clusterProfiler)
library(org.Gg.eg.db)

# 0.1 converting into bsseq object
BSobj = makeBSseqData( DF, c(names(DF)))

# 0.2 creation of the design matrix
design <- data.table()
for (i in 1:length(DF)){
     nom <- names(DF[i])
     cage <- substr(nom, nchar(nom), nchar(nom))
     w <- substr(nom, nchar(nom)-4, nchar(nom)-2)
     pair <- sub("^[A-Za-z](.*?)-.*", "\\1", nom)  
     design <- rbind(design, data.table(cage = cage, w = w, pair = pair))
}

# 1. création DML et DMR
DMLfit_interaction = DMLfit.multiFactor(BSobj, design=design, formula=~cage+w+pair+cage:w)
DML_WxCAGE = DMLtest.multiFactor(DMLfit_interaction, coef=4)
DMR_WxCAGE = callDMR(DMLtest_WxCAGE, p.threshold=0.05)


# 2. Transformer DML en objet GRanges
dml_hypo <- DML_WxCAGE %>% filter(stat < 0)
dml_hyper <- DML_WxCAGE %>% filter(stat > 0)

gr_dml_hypo <- GRanges(
  seqnames = dml_hypo$chr,
  ranges = IRanges(start = dml_hypo$pos, end = dml_hypo$pos),
  pval = dml_hypo$fdrs
)
gr_dml_hyper <- GRanges(
  seqnames = dml_hyper$chr,
  ranges = IRanges(start = dml_hyper$pos, end = dml_hyper$pos),
  pval = dml_hyper$fdrs
)

# 3. Importer le GFF
gff <- import.gff3("/work/project/geronimo/WP1/Jonathan/useful/annotation_promoteur_fixed.gff")
type_list <- c("gene", "exon", "five_prime_UTR", "three_prime_UTR", "promoter", "repeat", "TSS")
gff_filtered <- gff[mcols(gff)$type %in% type_list]

# 4. Annoter : trouver les overlaps
hits_hypo <- findOverlaps(gr_dml_hypo, gff_filtered)
hits_hyper <- findOverlaps(gr_dml_hyper, gff_filtered)

# 5.1 Joindre les annotations
annotated_dml_hypo <- as_tibble(gr_dml_hypo[queryHits(hits_hypo)])
annotated_dml_hypo$feature_type <- gff_filtered$type[subjectHits(hits_hypo)]
annotated_dml_hypo$gene_id <- gff_filtered$gene_id[subjectHits(hits_hypo)] 
annotated_dml_hyper <- as_tibble(gr_dml_hyper[queryHits(hits_hyper)])
annotated_dml_hyper$feature_type <- gff_filtered$type[subjectHits(hits_hyper)]
annotated_dml_hyper$gene_id <- gff_filtered$gene_id[subjectHits(hits_hyper)]


# 5.2 séléction des annot
annotated_dml_hypo_promoteur <- annotated_dml_hypo %>% filter(feature_type == "promoter") %>% filter(pval < 0.05)
annotated_dml_hypo_gene <- annotated_dml_hypo %>% filter(feature_type == "gene") %>% filter(pval < 0.05)
annotated_dml_hyper_promoteur <- annotated_dml_hyper %>% filter(feature_type == "promoter") %>% filter(pval < 0.05)
annotated_dml_hyper_gene <- annotated_dml_hyper %>% filter(feature_type == "gene") %>% filter(pval < 0.05)


# 6. Analyse GO (Biological Process)

ego_hyper_p <- enrichGO(
  gene = annotated_dml_hyper_promoteur$gene_id,
  OrgDb = org.Gg.eg.db,
  keyType = "ENSEMBL", 
  ont = "BP",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

ego_hypo_p <- enrichGO(
  gene = annotated_dml_hypo_promoteur$gene_id,
  OrgDb = org.Gg.eg.db,
  keyType = "ENSEMBL", 
  ont = "BP",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

ego_hyper_g <- enrichGO(
  gene = annotated_dml_hyper_gene$gene_id,
  OrgDb = org.Gg.eg.db,
  keyType = "ENSEMBL", 
  ont = "BP",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

ego_hypo_g <- enrichGO(
  gene = annotated_dml_hypo_gene$gene_id,
  OrgDb = org.Gg.eg.db,
  keyType = "ENSEMBL", 
  ont = "BP",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)

# 7. Visualiser les résultats d'enrichissement

# HYPER P
df <- as.data.frame(ego_hyper_p)

df$GeneRatioNumeric <- sapply(df$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})
df <- df[order(-df$GeneRatioNumeric), ]


dotplot_ego_hyper_p <- ggplot(df[1:10, ], aes(x = GeneRatioNumeric, y = reorder(Description, Count), size = Count)) +
  geom_point(color = "black") +
  theme_minimal() +
  labs(title = "Enrichissement hyperméthylation promoteurs",
       x = "Gene ratio",
       y = "GO Term") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/dotplot_GO_terms_hyper_p.png", dotplot_ego_hyper_p)

# HYPO P 

df <- as.data.frame(ego_hypo_p)

df$GeneRatioNumeric <- sapply(df$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})
df <- df[order(-df$GeneRatioNumeric), ]


dotplot_ego_hypo_p <- ggplot(df[1:10, ], aes(x = GeneRatioNumeric, y = reorder(Description, Count), size = Count)) +
  geom_point(color = "black") +
  theme_minimal() +
  labs(title = "Enrichissement hypométhylation promoteurs",
       x = "Gene ratio",
       y = "GO Term") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/dotplot_GO_terms_hypo_p.png", dotplot_ego_hypo_p)

# HYPER G

df <- as.data.frame(ego_hyper_g)

df$GeneRatioNumeric <- sapply(df$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})
df <- df[order(-df$GeneRatioNumeric), ]


dotplot_ego_hyper_g <- ggplot(df[1:10, ], aes(x = GeneRatioNumeric, y = reorder(Description, Count), size = Count)) +
  geom_point(color = "black") +
  theme_minimal() +
  labs(title = "Enrichissement hyperméthylation gène",
       x = "Gene ratio",
       y = "GO Term") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/dotplot_GO_terms_hyper_g.png", dotplot_ego_hyper_g)

# HYPO G

df <- as.data.frame(ego_hypo_g)

df$GeneRatioNumeric <- sapply(df$GeneRatio, function(x) {
  parts <- strsplit(x, "/")[[1]]
  as.numeric(parts[1]) / as.numeric(parts[2])
})
df <- df[order(-df$GeneRatioNumeric), ]


dotplot_ego_hypo_g <- ggplot(df[1:10, ], aes(x = GeneRatioNumeric, y = reorder(Description, Count), size = Count)) +
  geom_point(color = "black") +
  theme_minimal() +
  labs(title = "Enrichissement hypométhylation gène",
       x = "Gene ratio",
       y = "GO Term") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/dotplot_GO_terms_hypo_g.png", dotplot_ego_hypo_g)

# 8.network visualisation GO

network_plot <- cnetplot(ego_hypo_g, showCategory = 10, layout = "fr", colorEdge = TRUE, categorySize = "GeneRatio")
ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/network_GO_hypo_g.png", plot = network_plot)

network_plot <- cnetplot(ego_hyper_p, showCategory = 10, layout = "fr", colorEdge = TRUE, categorySize = "GeneRatio")
ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/network_GO_hyper_p.png", plot = network_plot)

network_plot <- cnetplot(ego_hypo_p, showCategory = 10, layout = "fr", colorEdge = TRUE, categorySize = "GeneRatio")
ggsave("/work/project/geronimo/WP1/Jonathan/final/plot/network_GO_hypo_p.png", plot = network_plot)

# 8. voix de signalisation KEGG

dml_hypo_entrez <- bitr(annotated_dml_hypo$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Gg.eg.db)
dml_hyper_entrez <- bitr(annotated_dml_hyper$gene_id, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Gg.eg.db)
