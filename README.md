# Differential analysis with DSS
<p align="center"> 
    <img src="logo/UT3.jpg" width="300" style="margin: 30px;">
    <img src="logo/INRAE.png" width="300" style="margin: 30px;">
    <img src="logo/geronimo.png" width="300" style="margin: 30px;">
    <img src="logo/genphyse.png" width="300" style="margin: 30px;">
</p>

## 📑 Table of Contents

- [1. Environment Setup](#1-setup-your-environnent)
- [2. Pipeline for Data Processing](#2-pipeline-for-data-processing)
- [3. Data Formatting and Statistics](#3-data-formating-and-statistics)
- [4. Differential Analysis](#4-differential-analysis)

## 1. Setup your environnent
**first download docker**
  
*if you use linux*
```sh
#only if you don t have already curl
sudo apt-get update
sudo apt-get install curl
#get docker
curl -fsSL https://get.docker.com/ | sh
```
*if you use windows click on the link [get docker](https://desktop.docker.com/win/main/amd64/Docker%20Desktop%20Installer.exe?utm_source=docker&utm_medium=webreferral&utm_campaign=docs-driven-download-win-amd64)*

**dowload the docker image to obtain all software and packages for analysis (our you can creat it by your own by clicking on [tuto env](docker/README.md)**
```sh
docker pull jchalaye/dss:analysis
```
**run the container and change the path to your working directory**
```sh
docker run -it --rm -v "path/to/your/working/directory:/work" dss:analysis
```
## 2. Pipeline for data processing
**Go in your work directory**
```sh
cd path/to/your/working/directory
```
**Copy the pipeline**
```sh
git clone https://github.com/seynard/pipeline_RRBS.git
```
**Go into the git repository:**
```sh
cd pipeline_RRBS
```

**Make sure you have the last version**
```sh
git pull
```
**Go into the nextflow directory**
```sh
cd nextflow_pipeline
```
**Edit the run_pipeline.sh file**
  
**Usage:**

<div style="border:1px solid #ccc; padding:10px; border-radius:5px;">

   - The first lines starting with ```SBATCH``` are slurm parameters
       
       - ```SBATCH -J``` is used to set your job name
    
       - ```SBATCH -o``` is used to set the output/log file
    
       - ```SBATCH -e``` is used to set the error file
    
       - ```SBATCH --mem=2G``` is used to set the memory for the main job (the pipeline will automatically set memory and cpus parameters for the jobs it launches)

   - The following lines with ```modules``` are used to load necessary tools, here only nextflow is necessary
  
   - The next lines corresponds to the nextflow command line to lauch the pipeline. The first lines is common nextflow option (for all pipelines) while the following ones are specific parameters for this pipeline. **Warning:** all this lines make one command line and thus **must finish with ```\``` except the last one** 

     - Common nextflow options:
     
       - ```-profile cluster``` is used to set execution environement. Use ```cluster``` value for genotoul execution, this will allow to automatically scale memory and cpus requested as well as define modules used during execution.
      
       - ```-resume``` allows nextflow to restart from prevously finished process if an error occurred or process are changed. **Warning: this is based on nextflow work directory if it is removed this will not work**
      
       - ```-with-report nf_report.html``` outputs a html report with process command and information on memory/cpus usage.

      - Specific parameters:
       
       - ```--samplesheet "/path/to/samplesheet.tsv"``` the samplesheet is a tabulation separated value file (.tsv). It must contains the header line with the 4 following columns: sample, read1, read2, umi.
           
           - ```sample``` sample id
        
           - ```read1``` path to the first read in the pair fastq
        
           - ```read2``` path to the second read in the pair fastq
        
           -  ```umi``` path to the umi fastq
       
    
       - ```--refdir "path/to/refefrence directory"``` path to the reference **directory without the trailing ```/```**.
    
       - ```--reference "Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"``` name of the fasta reference file.
    
       - ```--outdir "results_test"``` path to the general output directory. **This directory will only contain links to output file** as well as other non long term output.
    
       - ```--outbam "/work/project/geronimo/data/azenta/WP1_chicken_blood_novo/bam"``` path to the kept bam directory, where output bam will be copied.
    
       - ```--outvcf "/work/project/geronimo/data/azenta/WP1_chicken_blood_novo/vcf"``` path to the kept vcf directory, where output vcf will be copied.
    
       - ```--outvcf_filtered "/work/project/geronimo/data/azenta/WP1_chicken_blood_novo/vcf_processed"``` path to the kept filtered vcf directory, where output filtered vcf will be copied.
    
       - ```--outbed "/work/project/geronimo/data/azenta/WP1_chicken_blood_novo/bed"``` path to the kept bed directory, where output bed will be copied.
    
       - ```--outbed_filtered "/work/project/geronimo/data/azenta/WP1_chicken_blood_novo/bed_processed"``` path to the kept filtered bed directory, where output filtered bed will be copied.
    
       - ```--outreport "/work/project/geronimo/data/azenta/WP1_chicken_blood_novo/report"``` path to the kept report directory, where output report and logs will be copied.
    
       - ```--min_SNP_qual 10``` numeric value used to filter SNP on quality. (default: 10)
    
       -  ```--min_SNP_depth 10``` numeric value used to filter SNP on minimum depth. (default: 10)
    
       -  ```--max_SNP_depth "quantile"``` mode `"quantile"` or numeric value `--max_SNP_depth 150` used to filter SNP on maximum depth. (default: "quantile").

      If `"quantile"` mode is used depth greater than the 99.9th quantile will be filtered out. The quantile will be computed for each SAMPLE INDIVIDUALLY.
      If a numeric value is used depth greater than the given value will be filtered out for ALL SAMPLES. 
    
       -  ```--min_CG_depth 10``` numeric value used to filter CpG on minimum depth. (default: 10)
    
       -  ```--max_CG_depth "quantile"``` mode `"quantile"` or numeric value `--max_CG_depth 150` used to filter CpG on maximum depth. (default: "quantile").
         
      If `"quantile"` mode is used, depth greater than the 99.9th quantile will be filtered out. The quantile will be computed for each SAMPLE INDIVIDUALLY.
      If a numeric value is used depth greater than the given value will be filtered out for ALL SAMPLES.
</div>
  

**Lauch the pipeline**
```sh
sbatch run_pipeline.sh
```

## 3. Data formating and statistics

**go to working directory**
```sh
cd /work/project/geronimo/WP1/Jonathan/final
```

**creat directory of interest**
```sh
mkdir data depth variant DSS
```

**copy all scripts to your working directory**
```sh
git clone https://github.com/Jonathano3/DSSdiffanalysis
```

**run ["stat_unique_CpG.sh"](script/data_formating/stat_unique_CpG.sh) to collect statistics on unique CpGs during the process**
```sh
sbatch /work/project/geronimo/WP1/Jonathan/final/DSSdiffanalysis/script/data_formating/stat_unique_CpG.sh
```
*visualisation of the running script : ***stat_unique_CpG.sh***.*
```sh
dir_in=/work/project/geronimo/data/azenta/WP1_chicken_blood_novo
dir_out=/work/project/geronimo/WP1/Jonathan/final
fasta=/work/project/geronimo/data/azenta/fasta/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa
module load devel/python/Python-3.9.18
module load bioinfo/deepTools/3.5.2
module load bioinfo/samtools/1.20
module load bioinfo/bedtools/2.30.0
module load devel/Miniconda/Miniconda3
module load bioinfo/BISCUIT/1.4.0

# count CpG unique before filters (1. merge strands, 2. check for unique positions)
#1
mkdir -p ${dir_out}/bed_merged ${dir_out}/bed_processed ${dir_out}/count
cp ${dir_in}/bed_processed/* ${dir_out}/bed_processed

lf=($(ls ${dir_in}/bed/*.bed.gz))
for i in ${lf[@]}
do
name=$(echo ${i} | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
echo $i
sbatch --mem=20G --wrap="biscuit mergecg ${fasta} ${i} > ${dir_out}/bed_merged/${name}_merged.bed; gzip ${dir_out}/bed_merged/${name}_merged.bed"
done

#2
sbatch --wrap="python ${dir_out}/script/count_unique_CpG.py ${dir_out}/bed_merged  > ${dir_out}/count/count_merged.txt"

# check for unique positions after filter on min depth 10
sbatch --wrap="python ${dir_out}/script/count_unique_CpG.py ${dir_in}/bed_processed > ${dir_out}/count/count_min.txt"

# check for unique positions after removed low CpG in samples
sbatch --wrap="python ${dir_out}/script/count_unique_CpG.py ${dir_out}/bed_processed > ${dir_out}/count/count_lcpg"
```

**run the filter for a low number of CpGs per sample by executing this script ["filter_low_CpG.sh"](script/data_formating/filter_low_CpG.sh)**
```sh
sbatch /work/project/geronimo/WP1/Jonathan/final/DSSdiffanalysis/script/data_formating/filter_low_CpG.sh
```
*visualisation of the running script : ***filter_low_CpG.sh***.*
```sh
# removed bed_processed with low CpG count
seuil=200000
mkdir ${dir_out}/bed_processed/samples_low_CpG ${dir_out}/bed_merged/samples_low_CpG


for i in ${dir_out}/bed_processed/*.bed.gz
do
line_count=$(zcat "$i" | wc -l)
if [ "$line_count" -lt ${seuil} ];
then
nom=$(basename "$i" .bed.gz)
mv "${dir_out}/bed_processed/$nom"* "${dir_out}/bed_processed/samples_low_CpG/"
fi
done

#remove bed_merged with low CpG count
for file in ${dir_out}/bed_processed/samples_low_CpG/*bed.gz
do
name=$(basename ${file} | cut -d'.' -f1)
mv ${dir_out}/bed_merged/${name}* ${dir_out}/bed_merged/samples_low_CpG/
done
```
**you can check the correct structure of the data**
```sh
# verification merged files = 0 CpGs with 1 gap
compteur=0
for file in ${dir_out}/bed_merged/*.bed.gz
do
count=$(zcat "$file" | awk '{if ($3 - $2 == 1) print $0}' | wc -l)
compteur=$((compteur + count))
done
echo "Nombre total de lignes avec un ecart de 1 : $compteur"
```
**run this script ["variants_list.sh"](script/data_formating/variants_list.sh) to obtain a list of variants collected by BISCUIT**
```sh
sbatch /work/project/geronimo/WP1/Jonathan/final/DSSdiffanalysis/script/data_formating/variants_list.sh
```
*visualisation of the running script : ***variants_list.sh***.*
```sh
module load bioinfo/Bcftools/1.21
module load bioinfo/PLINK/2.00a4
module load bioinfo/VCFtools/0.1.16
dir_in=/work/project/geronimo/WP1/Jonathan/final
dir_out=/work/project/geronimo/WP1/Jonathan/final/variant
final=/work/project/geronimo/WP1/Jonathan/final/data

#extraction of samples of interest
samples=($(ls ${dir_in}/vcf_processed/*.vcf.gz))
echo -n > ${dir_out}/list_files.txt
for i in ${samples[@]}
do
nline=($(zcat ${i} | wc -l))
echo $nline >> ${dir_out}/nb_line_vcf.txt
echo $nline
if [ ${nline} -gt 120000 ] # here fix the minimum number of variants found in the raw vcf (if too little the individual might cause to discard many variants due to missing rate too high)
then 
echo ${i} >> ${dir_out}/list_files.txt
fi
done

#get statisctics
samples=($(ls ${dir_in}/vcf_processed/*.vcf.gz))
echo -n > ${dir_out}/nb_line_vcf.txt
for i in ${samples[@]}
do
nline=($(zcat ${i} | wc -l))
echo $nline >> ${dir_out}/nb_line_vcf.txt
done

#merge vcf
ulimit -n 2000; bcftools merge --file-list ${dir_out}/list_files.txt -Oz -o ${dir_out}/merged.vcf.gz
bcftools annotate --set-id +'%CHROM:%POS' ${dir_out}/merged.vcf.gz -Oz -o ${dir_out}/merged_rs.vcf.gz
rm ${dir_out}/merged.vcf.gz

#extraction of variant
cd ${dir_out}
plink2 --vcf merged_rs.vcf.gz --missing --allow-extra-chr --chr-set 39
plink2 --vcf merged_rs.vcf.gz --freq --allow-extra-chr --chr-set 39
MAF=0.05
MISS=0.9
QUAL=20
MIN_DEPTH=10
MAX_ALL=2
vcftools --gzvcf merged_rs.vcf.gz --remove-indels --max-alleles $MAX_ALL --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --minDP $MIN_DEPTH --recode --recode-INFO-all --stdout | bgzip -c > merged_filter.vcf.gz
bcftools stats merged_filter.vcf.gz > merged_filter.stats
bcftools query -H -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' merged_filter.vcf.gz > ${final}/snp_chicken_biscuit.txt
bzip2 ${final}/snp_chicken_biscuit.txt
```
**this script ["filter.r"](script/data_formating/filter.r) is used to create your final list of CpGs filtered with variants and missing values and obtain statistics**
```sh
sbatch /work/project/geronimo/WP1/Jonathan/final/DSSdiffanalysis/script/data_formating/run.sh
```
*visualisation of the running script : ***filter.r***.*
```sh
#load library
library(readxl)
library(data.table)
library(tidyverse)
library(reshape2)

#load metadata
data <- read_excel("/work/project/geronimo/data/azenta/metadata/WP1_blood_novo_1149hens.xlsx")

#get all file paths
lf=list.files(path='/work/project/geronimo/WP1/Jonathan/final/bed_processed/',pattern='bed.gz')
n<-grep('tbi',lf)
lf<-lf[-n]

# old filter on 1/1000 quantile per chromosome NO KEEP and name formatting with conditions (environment and period)

DF <- list()
dmax <- list()
for(i in lf){
  print(i)
  df <- fread(paste0('/work/project/geronimo/WP1/Jonathan/final/bed_processed/', i))
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
snp1<-fread('/work/project/geronimo/WP1/Jonathan/final/data/snp_chicken_biscuit.txt')
snp2<-fread('/work/project/geronimo/data/azenta/metadata/snp_chicken_novo.txt.gz')
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
fwrite(list_cpg,'/work/project/geronimo/WP1/Jonathan/final/data/list_cpg_final_sb2.txt',col.names=F,row.names=F,quote=F)

#save filter steps (depth)
fwrite(DF_depth,'/work/project/geronimo/WP1/Jonathan/final/depth/DF_depth.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(DF_sub_depth,'/work/project/geronimo/WP1/Jonathan/final/depth/DF_sub_depth.txt',col.names=T,row.names=F,quote=F, sep = "\t")
fwrite(DF_sub2_depth,'/work/project/geronimo/WP1/Jonathan/final/depth/DF_sub2_depth.txt',col.names=T,row.names=F,quote=F, sep = "\t")
```
**this script ["create_table.r"](script/data_formating/create_table.r) is used to create your final data table (methylation and depth) merged with all the CpGs found.**
```sh
R path/to/your/working/directory/DSSdiffanalysis/script/data_formating/create_table.r
```
*visualisation of the running script : ***create_table.r***.*
```sh
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
```
## 4. Differential analysis