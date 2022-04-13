FROM  continuumio/miniconda3:4.10.3
#FROM python:3.9-bullseye


RUN apt-get update
RUN apt-get install -y \
    autoconf \
    build-essential \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    wget \
    zlib1g-dev

WORKDIR /opt

# GeneFuse
RUN wget https://github.com/OpenGene/GeneFuse/archive/refs/tags/v0.6.1.tar.gz -O GeneFuse.tar.gz && \
    tar -xf GeneFuse.tar.gz && \
    rm GeneFuse.tar.gz

RUN cd GeneFuse-0.6.1 && \
    make && \
    cp genefuse /usr/local/bin/

# Lumpy
RUN conda create -n lumpy-sv \
    -c bioconda -c conda-forge -c default \
    lumpy-sv=0.3.1

# HmnFusion
COPY *whl /opt/
RUN pip install --no-cache-dir /opt/*whl

ENTRYPOINT ["hmnfusion"]
