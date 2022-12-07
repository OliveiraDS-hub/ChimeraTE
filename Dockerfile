FROM continuumio/miniconda

WORKDIR /chimeraTE

COPY . .

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda

RUN conda install -n base python=3.6.15 bedtools=2.30 blast=2.13 bowtie=1.3.1 cufflinks=2.2.1 curl express=1.5.1 hmmer=3.3.2 htslib=1.10.2 numpy=1.19 openjdk=11.0.1 perl=5.32.1 repeatmasker=4.1.2 salmon=1.7.0 samtools=1.10 seqtk=1.3 trimmomatic=0.39 trinity=2.11.0 wget=1.20.3 -y

RUN chmod -R 777 /chimeraTE *
