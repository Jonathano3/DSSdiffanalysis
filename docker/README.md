## Creat your own image docker

**go in your working directory**
```sh
cd path/to/your/working/directory
```
**create a file named “Dockerfile”**
```sh
touch Dockerfile
```
**edit it with the following code**
```sh
#source image
FROM continuumio/miniconda3

COPY environment.yaml /tmp/environment.yml

# Create a new conda environment with the specified software
RUN conda env create -f /tmp/environment.yml
# Use the new environment as the default
RUN echo "conda activate ENV_DSSdiffanalysis_data_formating" >> ~/.bashrc

# Use bash as the default shell
SHELL ["/bin/bash", "-c"]

# Install dependencies
RUN apt-get update && apt-get install -y \
    libfreetype6-dev \
    libfontconfig1-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN source activate ENV_DSSdiffanalysis_data_formating && \
    R -e "install.packages('BiocManager', ask = FALSE, repos='https://cloud.r-project.org/')" && \
    R -e "BiocManager::install(version = '3.20')" && \
    R -e "BiocManager::install(c('reshape2', 'data.table', 'tidyverse', 'readxl', 'ggplot2'))" && \
    R -e "BiocManager::install(c('GenomicRanges', 'AnnotationDbi', 'rtracklayer', 'biomaRt', 'DSS', 'clusterProfiler'))"

#set the working directory
WORKDIR /work

# set the entrypoint to bash
ENTRYPOINT ["/bin/bash"]
```

**create another file named “environment.yaml”**
```sh
touch environment.yaml
```
**edit it with the following code**
```sh
name: ENV_DSSdiffanalysis_data_formating
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.9.18
  - r-base=4.4.0
  - samtools=1.20
  - bedtools=2.30.0
  - deeptools=3.5.2
  - biscuit=1.4.0
  - Bcftools=1.21
  - PLINK=2.00a4
  - VCFtools=0.1.16
```
**creat your docker image**
```sh
docker build -t dss:analysis .
```
**you can run it and follow [the next step of the analysis (2)](../README.md#2-pipeline-for-data-processing)**
```sh
docker run -it --rm -v "path/to/your/working/directory:/work" dss:analysis
```