#!/bin/bash
#SBATCH --job-name=stat_unique_CpG


dir_in=path/to/your/working/directory/data
dir_out=path/to/your/working/directory
script=path/to/your/working/directory/DSSdiffanalysis/script
fasta=path/to/your/working/directory/data/fasta/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa


# count CpG unique before filters (1. merge strands, 2. check for unique positions)
#1
mkdir -p ${dir_out}/bed_merged ${dir_out}/bed_processed ${dir_out}/count
cp ${dir_in}/bed_processed/* ${dir_out}/bed_processed

lf=($(ls ${dir_in}/bed/*.bed.gz))
for i in ${lf[@]}
do
name=$(echo ${i} | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
echo $i
biscuit mergecg ${fasta} ${i} > ${dir_out}/bed_merged/${name}_merged.bed; gzip ${dir_out}/bed_merged/${name}_merged.bed
done

#2
python ${script}/count_unique_CpG.py ${dir_out}/bed_merged  > ${dir_out}/count/count_merged.txt

# check for unique positions after filter on min depth 10
python ${script}/count_unique_CpG.py ${dir_in}/bed_processed > ${dir_out}/count/count_min.txt

# check for unique positions after removed low CpG in samples
python ${script}/count_unique_CpG.py ${dir_out}/bed_processed > ${dir_out}/count/count_lcpg
