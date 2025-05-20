#!/bin/bash
#SBATCH --job-name=filter_low_CpG

# removed bed_processed with low CpG count
seuil=200000
dir_out=path/to/your/working/directory
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
