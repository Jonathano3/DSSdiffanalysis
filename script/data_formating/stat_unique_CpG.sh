#!/bin/bash
#SBATCH --job-name=stat_unique_CpG


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
