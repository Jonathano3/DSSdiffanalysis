#load library
library(readxl)
library(data.table)
library(tidyverse)
library(reshape2)

#load metadata
data <- read_excel("path/to/your/working/directory/data/WP1_blood_novo_1149hens.xlsx")

#get all file paths
lf=list.files(path='path/to/your/working/directory/bed_processed/',pattern='bed.gz')
n<-grep('tbi',lf)
lf<-lf[-n]

# old filter on 1/1000 quantile per chromosome NO KEEP and name formatting with conditions (environment and period)

DF <- list()
dmax <- list()
for(i in lf){
  print(i)
  df <- fread(paste0('path/to/your/working/directory/bed_processed/', i))
  colnames(df) <- c('chr', 'start', 'stop', 'meth', 'depth', 'info')
  df$rs <- paste0(df$chr, ':', df$start, '_', df$stop)

  #dmax[[i]] <- list()
  #for (chrom in unique(df$chr)) {
  #  df_chr <- df[df$chr == chrom, ]
  #  dmax[[i]][[chrom]] <- quantile(df_chr$depth, probs = 0.999)  
  #}
  df_filtered <- df
  #for (chrom in unique(df$chr)) {
  #  seuil <- dmax[[i]][[chrom]]
  #  df_filtered <- df_filtered[!(df_filtered$chr == chrom & df_filtered$depth >= seuil), ]
  #}
  name = sub("\\.filtered\\.bed\\.gz$", "", i)
  w = paste0("-",trimws(data$ageFactor[data$sampleIdRRBS == name]))
  cage = paste0("-",trimws(data$isCageUntil55w[data$sampleIdRRBS == name]))
  DF[[i]] <- df_filtered
  DF[[i]]$ind <- paste0(name,w,cage)
}

DF<-do.call(rbind,DF)

##table creation
# subdivide the data frame to avoid the problem
individus_uniques <- unique(DF$ind)
pre_moitier_individus <- tail(unique(DF$ind), 500)
autre_moitie_individus <- individus_uniques[!individus_uniques %in% pre_moitier_individus]

df_filtered_depth <- DF[ind %in% pre_moitier_individus]
df_filtered2_depth <- DF[ind %in% autre_moitie_individus]
DF_depth_1 <- spread(df_filtered_depth[,c('rs','depth','ind')],ind,depth)
DF_depth_2 <- spread(df_filtered2_depth[,c('rs','depth','ind')],ind,depth)
DF_depth<- merge(DF_depth_1, DF_depth_2, by = "rs", all = TRUE)


### Variant filter
#load varaint files
snp1<-fread('path/to/your/working/directory/data/snp_chicken_biscuit.txt')
snp2<-fread('path/to/your/working/directory/data/snp_chicken_novo.txt.gz')
colnames(snp1)<-c('chr','pos','rs','ref','alt')
colnames(snp2)<-c('chr','pos','rs','ref','alt')
snp<-rbind(snp1,snp2)
snp$rs<-NULL

#select snps of interest
snp<-unique(snp)
snp<-subset(snp,snp$ref%in%c('C','G'))

#fitting with correct start format
a<-colsplit(DF_depth$rs,':',c('chr','pos'))
DF_depth$chr<-a$chr
b<-colsplit(a$pos,'_',c('start','end'))
DF_depth$pos<-b$start+1
chrom<-unique(snp$chr)

#filtering variants
DF_sub_depth<-list()
for (i in chrom) {
  x <- DF_depth[DF_depth$chr == i, ]
  snp_i <- snp[snp$chr == i, ]
  snp_C <- snp_i$pos[snp_i$ref == "C"]
  snp_G <- snp_i$pos[snp_i$ref == "G"]
  y <- subset(x, !(x$pos%in%c(snp_C, snp_C + 1,snp_G , snp_G - 1)))
  DF_sub_depth[[i]] <- y
}
DF_sub_depth<-do.call(rbind,DF_sub_depth)

###missing filter
DF_sub_depth$miss<-rowSums(is.na(DF_sub_depth))/(ncol(DF_sub_meth)-2)
DF_sub2_meth<-DF_sub_meth[DF_sub_meth$miss<0.25,]

#statistics depth
n<-colnames(DF_sub2_depth)[!colnames(DF_sub2_depth)%in%c('chr','pos','rs','miss')]
DF_sub2_long_depth<-gather(DF_sub2_depth,ind,depth,n)

n1 <- nrow(DF_depth)
n2 <- nrow(DF_sub_depth)
n3 <- nrow(DF_sub2_depth)
moyenne <- mean(DF_sub2_long_depth$depth, na.rm = TRUE)
statistiques <- c(
  paste("Nombre de lignes dans DF_depth :", n1),
  paste("Nombre de lignes dans DF_sub_depth :", n2),
  paste("Nombre de lignes dans DF_sub2_depth :", n3),
  paste("Profondeur moyenne dans DF_sub2_long_depth :", moyenne)
)

#write statistics
writeLines(statistiques, "/work/project/geronimo/WP1/Jonathan/final/data/statistiques_depth.txt")

#write list CpGs of interest
list_cpg<-DF_sub2_depth[,c('chr','pos')]
list_cpg$end<-DF_sub2_depth$pos+1
fwrite(list_cpg,'path/to/your/working/directory/data/list_cpg_final_sb2.txt',col.names=F,row.names=F,quote=F)

#save filter steps (depth)
fwrite(DF_depth,'path/to/your/working/directory/depth/DF_depth.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(DF_sub_depth,'path/to/your/working/directory/depth/DF_sub_depth.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(DF_sub2_depth,'path/to/your/working/directory/depth/DF_sub2_depth.txt',col.names=T,row.names=F,quote=F, sep = "\t")
