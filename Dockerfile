FROM ubuntu:20.04

LABEL maintainer="Daniel S. Oliveira <daniel.sdo2015@gmail.com>"
LABEL version="ChimeraTE v1.3"
LABEL date="May 1st 2025"


ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C

WORKDIR /opt/chimerate

RUN apt update && apt install -y \
    pkg-config \
    libhdf5-dev \
    libhdf5-serial-dev \
    libhdf5-103 \
    ca-certificates \
    apt-utils \
    autoconf \
    automake \
    cmake \
    gcc \
    wget \
    build-essential \
    software-properties-common \
    tar \
    unzip \
    zlib1g-dev \
    sudo \
    git-core \
    locales \
    python3-dev \
    python3-pip \
    libbz2-dev \
    liblzma-dev \
    libncurses-dev \
    libfile-which-perl \
    libjson-perl \
    liburi-perl \
    liblwp-useragent-determined-perl \
    libtext-soundex-perl \
    perl \
    default-jdk
    
RUN pip3 install --upgrade pip && \
    python3 -m pip install --no-cache-dir \
    biopython termcolor==1.1.0 pandas==1.1.5 python-dateutil==2.8.2 setuptools==59.6.0 \
    Cython==0.29.16 --no-binary=h5py h5py numpy==1.19.5 pyranges==0.1.4

# STAR
RUN wget https://github.com/alexdobin/STAR/releases/download/2.7.10b/STAR_2.7.10b.zip && \
    unzip STAR_2.7.10b.zip -d /usr/local && \
    rm STAR_2.7.10b.zip

# Bowtie2
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/bowtie2-2.5.1-linux-x86_64.zip && \
    unzip bowtie2-2.5.1-linux-x86_64.zip -d /usr/local && \
    rm bowtie2-2.5.1-linux-x86_64.zip

# Cufflinks
RUN wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz && \
    tar -xzf cufflinks-2.2.1.Linux_x86_64.tar.gz -C /usr/local && \
    rm cufflinks-2.2.1.Linux_x86_64.tar.gz

# Express
RUN wget https://pachterlab.github.io/eXpress/downloads/express-1.5.1/express-1.5.1-linux_x86_64.tgz && \
    tar -xzf express-1.5.1-linux_x86_64.tgz -C /usr/local && \
    rm express-1.5.1-linux_x86_64.tgz

# Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz && \
    tar -xzf bedtools-2.30.0.tar.gz && cd bedtools2 && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    make && make install && \
    cd .. && rm -rf bedtools2*

# Seqtk
RUN wget https://github.com/lh3/seqtk/archive/refs/tags/v1.3.zip && \
    unzip v1.3.zip && cd seqtk-1.3 && make && cp seqtk /usr/local/bin && \
    cd .. && rm -rf seqtk-1.3 v1.3.zip

# Samtools
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.7/samtools-1.7.tar.bz2 && \
    tar -xf samtools-1.7.tar.bz2 && cd samtools-1.7 && ./configure && make && make install && \
    cd .. && rm -rf samtools-1.7*

# Jellyfish
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.3.0/jellyfish-linux -O /usr/local/bin/jellyfish && \
    chmod +x /usr/local/bin/jellyfish

# Salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz && \
    tar -xzf salmon-1.10.0_linux_x86_64.tar.gz && \
    mv salmon-latest_linux_x86_64 /usr/local/salmon && \
    rm salmon-1.10.0_linux_x86_64.tar.gz

# Trinity
RUN wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/v2.9.1/trinityrnaseq-v2.9.1.FULL.tar.gz && \
    tar -xzf trinityrnaseq-v2.9.1.FULL.tar.gz && \
    cd trinityrnaseq-v2.9.1 && make && \
    cd .. && mv trinityrnaseq-v2.9.1 /usr/local/trinity && \
    rm trinityrnaseq-v2.9.1.FULL.tar.gz

# BLAST
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.14.1/ncbi-blast-2.14.1+-x64-linux.tar.gz && \
    tar -xzf ncbi-blast-2.14.1+-x64-linux.tar.gz -C /usr/local && \
    rm ncbi-blast-2.14.1+-x64-linux.tar.gz

# RepeatMasker dependencies
RUN wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64 -O /usr/local/bin/trf && chmod +x /usr/local/bin/trf

RUN wget http://eddylab.org/software/hmmer/hmmer-3.3.2.tar.gz && \
    tar -xzf hmmer-3.3.2.tar.gz && cd hmmer-3.3.2 && ./configure && make && make install && \
    cd .. && rm -rf hmmer-3.3.2*

RUN wget https://www.repeatmasker.org/rmblast/rmblast-2.13.0+-x64-linux.tar.gz && \
    tar -xzf rmblast-2.13.0+-x64-linux.tar.gz -C /usr/local && \
    rm rmblast-2.13.0+-x64-linux.tar.gz

RUN wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.2-p1.tar.gz && \
    tar -xzf RepeatMasker-4.1.2-p1.tar.gz && \
    mv RepeatMasker /usr/local/RepeatMasker && \
    rm RepeatMasker-4.1.2-p1.tar.gz && \
    /usr/local/RepeatMasker/configure -default_search_engine=rmblast \
        -rmblast_dir=/usr/local/rmblast-2.13.0/bin \
        -trf_prgm=/usr/local/bin/trf

# Add all relevant tools to PATH
ENV PATH="$PATH:/usr/local/STAR_2.7.10b/Linux_x86_64:/usr/local/bowtie2-2.5.1-linux-x86_64:/usr/local/cufflinks-2.2.1.Linux_x86_64:/usr/local/trinity:/usr/local/express-1.5.1-linux_x86_64:/usr/local/ncbi-blast-2.14.1+/bin:/usr/local/salmon/bin:/usr/local/RepeatMasker"

# Copy your script into the container
COPY projects/ ./projects/
COPY util/ ./util/
COPY scripts/ ./scripts/
COPY chimTE_mode1.py .
COPY chimTE_mode2.py .

# Ensure all Python files are executable (optional)
RUN chmod +x chimTE_mode*.py

# Add all relevant tools to PATH (same as before)
ENV PATH="$PATH:/usr/local/STAR_2.7.10b/Linux_x86_64:/usr/local/bowtie2-2.5.1-linux-x86_64:/usr/local/cufflinks-2.2.1.Linux_x86_64:/usr/local/trinity:/usr/local/express-1.5.1-linux_x86_64:/usr/local/ncbi-blast-2.14.1+/bin:/usr/local/salmon/bin:/usr/local/RepeatMasker"

# Default command — defer to argument passed at runtime
ENTRYPOINT ["python3", "-u"]
CMD ["chimTE_mode1.py"]
