## Creat your own image docker

**go in your working directory**
```sh
cd path/to/your/working/directory
```
**create a file named “Dockerfile”**
```sh
touch Dockerfile
```
**edit it with the following code or download the file (Dockerfile)[../script/Dockerfile]**
```sh
#source image
FROM continuumio/miniconda3

# Copy the environment.yaml file and create the conda environment
COPY environment.yaml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
RUN echo "conda activate ENV_DSSdiffanalysis_data_formating" >> ~/.bashrc

# Activate the Conda environment each time a command is run
SHELL ["conda", "run", "-n", "ENV_DSSdiffanalysis_data_formating", "/bin/bash", "-c"]

# Install system dependencies required for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglu1-mesa-dev \
    libxt-dev \
    zlib1g-dev && \
    apt-get clean

# Install required R packages
RUN R -e "install.packages('BiocManager', ask = FALSE, repos='https://cloud.r-project.org/')" && \
    R -e "BiocManager::install(version = '3.20')" && \
    R -e "BiocManager::install(c('GenomicRanges', 'AnnotationDbi', 'rtracklayer', 'biomaRt', 'DSS', 'clusterProfiler'))" && \
    R -e "install.packages('data.table', version = '1.17.0', repos = c(CRAN='https://cloud.r-project.org/'))" && \
    R -e "install.packages('ggplot2', version = '3.5.2', repos = c(CRAN='https://cloud.r-project.org/'))" && \
    R -e "install.packages('tidyverse', version = '2.0.0', repos = c(CRAN='https://cloud.r-project.org/'))" && \
    R -e "install.packages('readxl', version = '1.4.5', repos = c(CRAN='https://cloud.r-project.org/'))" && \
    R -e "install.packages('reshape2', version = '1.4.4', repos = c(CRAN='https://cloud.r-project.org/'))"


# Set the working directory
WORKDIR /work

# Set the default entry point to bash shell
ENTRYPOINT ["/bin/bash"]
```

**create another file named “environment.yaml”**
```sh
touch environment.yaml
```
**edit it with the following code**
```sh
name: ENV_DSSdiffanalysis_data_formating
channel_priority: strict
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.9.18
  - r-base=4.4.0
  - samtools=1.20
  - bedtools=2.30.0
  - deeptools=3.5.2
  - biscuit=1.4.0
  - bcftools=1.21
  - plink2=2.0.0a.6.9
  - vcftools=0.1.16
```
**creat your docker image**
```sh
docker build -t dss:analysis .
```
**you can run it and follow [the next step of the analysis (2)](../README.md#2-pipeline-for-data-processing)**
```sh
docker run -it --rm -v "path/to/your/working/directory:/work" dss:analysis
```