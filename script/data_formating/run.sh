#!/bin/bash
#SBATCH --job-name=filter_low_CpG
#SBATCH --mem=300G

module load statistics/R/4.4.0
Rscript /work/project/geronimo/WP1/Jonathan/final/DSSdiffanalysis/script/data_formating/filter.r
