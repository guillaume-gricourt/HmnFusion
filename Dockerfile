FROM python:3.9-bullseye


RUN apt-get update
RUN apt-get install -y \
    autoconf \
    build-essential \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \
    libncurses5-dev \
    wget

WORKDIR /opt

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar -xf samtools-1.15.tar.bz2 && \
    cd samtools-1.15 && \
    ./configure && \
    make && \
    make install && \
    cd /opt \
    rm samtools-1.15.tar.bz2

# GeneFuse
RUN wget https://github.com/OpenGene/GeneFuse/archive/refs/tags/v0.8.0.tar.gz -O GeneFuse.tar.gz && \
    tar -xf GeneFuse.tar.gz && \
    rm GeneFuse.tar.gz

RUN cd GeneFuse-0.8.0 && \
    make && \
    make install && \
    cd /opt

# Lumpy
RUN wget https://github.com/arq5x/lumpy-sv/releases/download/0.3.0/lumpy-sv.tar.gz -O Lumpy.tar.gz && \
    tar -xf Lumpy.tar.gz && \
    rm Lumpy.tar.gz

RUN cd lumpy-sv && \
    make && \
    cp bin/* /usr/local/bin/ && \
    cp scripts/extractSplitReads_BwaMem /usr/local/bin/ && \
    cp scripts/lumpyexpress /usr/local/bin/ && \
    cd /opt && \
    rm -rf lumpy-sv

# HmnFusion
COPY *whl /opt/
RUN pip install --no-cache-dir /opt/*whl

ENTRYPOINT ["hmnfusion"]
