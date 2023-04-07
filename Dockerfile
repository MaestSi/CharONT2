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

RUN conda create -n CharONT_env bioconductor-biostrings r-kernsmooth 
RUN conda install -n CharONT_env emboss vsearch seqtk mafft minimap2 samtools=1.15 medaka trf NanoFilt bbmap racon
RUN conda create -n pycoQC_env python=3.6 pip

RUN /opt/conda/envs/pycoQC_env/bin/python -m pip install pycoQC
RUN ln -s /opt/conda/envs/pycoQC_env/bin/pycoQC /opt/conda/envs/CharONT_env/bin
WORKDIR /home/
