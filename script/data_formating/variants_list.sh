#!/bin/bash
#SBATCH --job-name=variant_list
dir_in=path/to/your/working/directory
dir_out=path/to/your/working/directory/variant
final=path/to/your/working/directory/data

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