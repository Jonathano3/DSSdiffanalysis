# Visualization
This directory will show you all the visualizations I made during my Master 2 internship.  
  
Warning :  
* All scripts and visualizations to be shared are exclusively **tailor-made for the WP1_chicken_blood_novo** project data and my path.
* All scripts are run locally

## Metadata analysis  
  
First, a visualization of the raw data distribution is performed with this R script : **[metadata.r](../script/visualization/metadata.r)**

<p align="center">
  <img src="../figure/metadata/boxplot_lecture_brute.png" width="32%" style="display: inline-block; margin: 0;">
  <img src="../figure/metadata/density_distribution.png" width="32%" style="display: inline-block; margin: 0;">
  <img src="../figure/metadata/read_barcode.png" width="32%" style="display: inline-block; margin: 0;">
</p>

## Benchmarking maximum threshold

The custom threshold is then visualized with the script : **[benchmarking.r](../script/visualization/metadata.r)**
<p align="center">
  <img src="../figure/benchmarking/plot_37.png" width="32%" style="display: inline-block; margin: 0;">
  <img src="../figure/benchmarking/plot_density_37.png" width="32%" style="display: inline-block; margin: 0;">
  <img src="../figure/benchmarking/densite_secondaire_chrom37.png" width="32%" style="display: inline-block; margin: 0;">
</p>  

<p align="center">
  <img src="../figure/benchmarking/ratio_CpG_taille.png" width="50%" style="display: inline-block; margin: 0;">
</p>

## View filters
Visualization of some statistics after the pipeline with the script : **[View_filter.r](../script/visualization/View_filter.r)**
<p align="center">
  <img src="../figure/filter/seuil_bed_processed.png" width="32%" style="display: inline-block; margin: 0;">
  <img src="../figure/filter/valeur_NA_camembert.png" width="32%" style="display: inline-block; margin: 0;">
</p>  
  
## Data characterization
Now that our data is clean, we can study it further with this script : **[characterization.r](../script/visualization/characterization.r)**

<p align="center">
  <img src="../figure/characterization/perc_features2_sb2.png" width="32%" style="display: inline-block; margin: 0;">
  <img src="../figure/characterization/perc_island_sb2.png" width="32%" style="display: inline-block; margin: 0;">
</p> 

<p align="center">
  <img src="../figure/characterization/S.D_delta_70w_90w.png" width="32%" style="display: inline-block; margin: 0;">
  <img src="../figure/characterization/S.D_delta_no_cage_cage.png" width="32%" style="display: inline-block; margin: 0;">
</p>

## DSS analysis