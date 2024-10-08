#------------------------------------------------------------------------------------#
#
#   Created by Alem Gusinac, last modified at 26-09-2024
# 
#   This Dockerfile creates a container suitable for r-base:latest and Qiime2
#   It builds on top of a slim ubuntu (70MB) base
# 
#   Current available setups by section:
#       1. qiime2
#       2. R: dada2
# 
#   Additional setups can be included after the latest section
#   It is important to include lib/pkgs for your setup in the first section
#
#------------------------------------------------------------------------------------#

FROM ubuntu:20.04

ARG USER=docker

# set environment without graphical interface
ENV DEBIAN_FRONTEND=noninteractive

#------------------------------------------------------------------------------------#
# 1. Installation of essential linux packages
#------------------------------------------------------------------------------------#

# Install dependencies
RUN apt-get update && apt-get install -y \
    wget curl \
    libtool \
    libcurl4-openssl-dev \ 
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    zlib1g-dev \
    libncurses5-dev \
    libgdbm-dev \
    libnss3-dev \
    libreadline-dev \
    libffi-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

#------------------------------------------------------------------------------------#
# 2. R-base:latest setup
#------------------------------------------------------------------------------------#
# Setting up R-base:latest, includes public key setup and addition of r-base version to apt depository
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" >> /etc/apt/sources.list
# apt update, upgrade and clean after above set-up is required!
# Otherwise it gives an error of broken pipes.
RUN apt-get update && apt-get upgrade -y && apt-get clean
RUN apt-get update && apt-get install -y r-base

# Copy requirements
COPY install2.r .
COPY installBioc.r .
COPY installGithub.r .

# Required package for install2.r
RUN R -e "install.packages('docopt', dependencies=TRUE)"

# R dependencies from CRAN
RUN Rscript install2.r --error --skipinstalled \
    BiocManager \
    optparse \
    remotes \
    data.table \
    foreach \
    parallel \
    doParallel \
    purrr \
    dplyr \
    && rm -rf /tmp/downloaded_packages

# R dependencies from BiocManager
RUN Rscript installBioc.r --error --skipinstalled \
    dada2 \
    biomformat \
    phyloseq \
    && rm -rf /tmp/downloaded_packages

#------------------------------------------------------------------------------------#
# Final steps
#------------------------------------------------------------------------------------#

RUN useradd -ms /bin/bash ${USER}
USER ${USER}