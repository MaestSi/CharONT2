FROM continuumio/miniconda3

########### set variables
ENV DEBIAN_FRONTEND noninteractive

########## generate working directories
RUN mkdir /home/tools

######### dependencies
RUN apt-get update -qq \
    && apt-get install -y \
    build-essential \
    wget \
    unzip \
    bzip2 \
    git \
    libidn11* \
    nano \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

############################################################ install CharONT
WORKDIR /home/tools/

RUN git clone https://github.com/MaestSi/CharONT.git
WORKDIR /home/tools/CharONT
RUN chmod 755 *

RUN sed -i 's/PIPELINE_DIR <- .*/PIPELINE_DIR <- \"\/home\/tools\/CharONT\/\"/' config_CharONT.R
RUN sed -i 's/MINICONDA_DIR <- .*/MINICONDA_DIR <- \"\/opt\/conda\/\"/' config_CharONT.R

RUN conda config --add channels r && \
conda config --add channels anaconda && \
conda config --add channels conda-forge && \
conda config --add channels bioconda

RUN conda create -n CharONT_env bioconductor-biostrings r-kernsmooth r-factoextra python=3.8
RUN conda install -n CharONT_env emboss
RUN conda install -n CharONT_env vsearch
RUN conda install -n CharONT_env seqtk
RUN conda install -n CharONT_env mafft
RUN conda install -n CharONT_env minimap2
RUN conda install -n CharONT_env trf
RUN conda install -n CharONT_env bbmap
RUN conda install -n CharONT_env nanofilt
RUN conda install -n CharONT_env samtools
RUN /opt/conda/envs/CharONT_env/bin/python -m pip install medaka

RUN conda create -n pycoQC_env python=3.8 pip racon
RUN /opt/conda/envs/pycoQC_env/bin/python -m pip install pycoQC
RUN ln -s /opt/conda/envs/pycoQC_env/bin/pycoQC /opt/conda/envs/CharONT_env/bin
RUN ln -s /opt/conda/envs/pycoQC_env/bin/racon /opt/conda/envs/CharONT_env/bin
WORKDIR /home/
