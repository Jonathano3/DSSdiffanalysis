library(readxl)
library(data.table)
library(tidyverse)
library(reshape2)

#load your metadata
data <- read_excel("/work/project/geronimo/data/azenta/metadata/WP1_blood_novo_1149hens.xlsx")

#collect path
lf=list.files(path='/work/project/geronimo/WP1/Jonathan/final/bed_merged/',pattern='bed.gz')
CpG <- fread("/work/project/geronimo/WP1/Jonathan/final/data/list_cpg_final_sb2.txt", sep=",")
colnames(CpG) <- c('chr', 'start', 'stop')

#merged files
DF <- list()
for (i in lf) {
    print(i)
    df <- fread(paste0('/work/project/geronimo/WP1/Jonathan/final/bed_merged/', i))


    colnames(df) <- c('chr', 'start', 'stop', 'meth', 'depth', 'info')
    
    name = sub("_merged\\.bed\\.gz$", "", i)
    w = paste0("-", trimws(data$ageFactor[data$sampleIdRRBS == name]))
    cage = paste0("-", trimws(data$isCageUntil55w[data$sampleIdRRBS == name]))

    df$start <- df$start + 1
    df$rs <- paste0(df$chr, ':', df$start, '_', df$stop)
    df_filtered <- merge(df, CpG, by = c("chr", "start", "stop"))
    
    DF[[i]] <- df_filtered
    DF[[i]]$ind <- paste0(name, w, cage)
}

DF<-do.call(rbind,DF)

#create merged files
DF_meth<-spread(DF[,c('rs','meth','ind')],ind,meth)
DF_depth<-spread(DF[,c('rs','depth','ind')],ind,depth)

#save files
fwrite(DF_meth,'/work/project/geronimo/WP1/Jonathan/final/data/all_meth_final_unique_CpG_sb2.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(DF_depth,'/work/project/geronimo/WP1/Jonathan/final/data/all_depth_final_unique_CpG_sb2.txt',col.names=T,row.names=F,quote=F, sep = "\t")