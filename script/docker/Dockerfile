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
