![image](https://github.com/OliveiraDS-hub/ChimeraTE/blob/main/image/logo.png)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54) <br />
[![DOI](https://zenodo.org/badge/470747321.svg)](https://zenodo.org/badge/latestdoi/470747321)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/1681)



---
# The pipeline
ChimeraTE is a pipeline to detect chimeric transcripts derived from genes and transposable elements (TEs). It has two running Modes:

- Mode 1 chimeric transcripts detection based upon exons and TE copies positions in the genome sequence; 

- Mode 2 chimeric transcripts detection regardless the genomic position, allowing the detection of chimeras from TEs that are not present in the referece genome, but with less sensitivity.


1. [Install](#installation)
   1. [Conda](#conda)
   2. [Singularity](#singularity)
   3. [Requirements](#requirements)
2. [Required data](#req_data)
3. [ChimeraTE Mode 1](#mode1)
    1. [Preparing your data](#prep_data)
    2. [Example data](#example_m1)
    3. [Output](#output_m1)
4. [ChimeraTE Mode 2](#mode2)
    1. [Preparing your data](#prep_mode2)
    2. [Example data](#example_m2)
    3. [Output](#output_m2)

---

## Install <a name="installation"></a>
### Conda <a name="conda"></a>
The installation may be easily done with conda. If you don't have conda installed in your machine, please follow [this tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

Once you have installed conda, you need to enable Bioconda channel with:
````conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
````

Then, all dependencies to run ChimeraTE can be easily installed in a new conda environment by using the chimeraTE.yml file:

Download repository from github:<br />````git clone https://github.com/OliveiraDS-hub/ChimeraTE.git````

Change to the ChimeraTE's folder:<br />````cd ChimeraTE````

Create chimeraTE environment with all dependencies:<br />````conda env create -f chimeraTE.yml````

Activate the new environment:<br />````conda activate chimeraTE````

Note: We advise you to return your condarc config to the default with:
```
conda config --remove channels bioconda
conda config --remove channels conda-forge
conda config --set channel_priority false
```
### Singularity <a name="singularity"></a>
Alternatively to conda, you can use [singularity](https://github.com/sylabs/singularity/releases/tag/v3.10.0) v3.10.0+ to build a container with all dependencies for ChimeraTE.

If you don't have ```sudo``` permissions:
```
singularity build --fakeroot chimeraTE.simg singularity.def
```
If you have ```sudo```:
```
sudo singularity build chimeraTE.simg singularity.def
```

Then, to run ChimeraTE:

```singularity exec chimeraTE.simg python3 chimTE_mode1.py --help```
<br />```singularity exec chimeraTE.simg python3 chimTE_mode2.py --help```<br />





### Requirements <a name="requirements"></a>

If you don't have conda, or you are working with an OS that some packages in the conda repo are not available, you can install them manually. It's important to highlight that all of them must be installed in your path.

- Python dependencies
  - [Python 3.6+](https://www.python.org/downloads/release/python-360/)
  - [Termcolor 1.1.0](https://pypi.org/project/termcolor/1.1.0/)
  - [Pybedtools 0.9.0](https://daler.github.io/pybedtools/main.html)
  - [Pandas 1.1.5](https://pandas.pydata.org/docs/getting_started/install.html)
  - [Dateutil 2.8.2](https://pypi.org/project/python-dateutil/)
  - [h5py 3.1.0](https://pypi.org/project/h5py/3.1.0/)
  - [numpy 1.19.5](https://numpy.org/install/)
  
- Softwares
  - [Bedtools 2.30.0](https://github.com/arq5x/bedtools2/releases/tag/v2.30.0)
  - [BLAST 2.13+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
  - [Bowtie 2.5.1](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/)
  - [Cufflinks 2.2.1](http://cole-trapnell-lab.github.io/cufflinks/install/)
  - [express 1.5.1](https://pachterlab.github.io/eXpress/overview.html)
  - [RepeatMasker 4.1.2](https://www.repeatmasker.org/RepeatMasker/)
  - [Trinity 2.9.1](https://github.com/trinityrnaseq/trinityrnaseq/releases/tag/v2.9.1)
  - [Seqtk 1.3](https://github.com/lh3/seqtk)
  - [Samtools 1.7](https://sourceforge.net/projects/samtools/files/samtools/1.7/)
  - [STAR 2.7.10](https://github.com/alexdobin/STAR)

## Required data <a name="req_data"></a>
In order to run ChimeraTE, the following files are required according to the running Mode: 

| Data | Mode 1 | Mode 2 | Mode 2 --assembly |
| -------- | -------- | -------- | -------- |
| Stranded paired-end RNA-seq - Fastq files     | X    |   X   |    X   |
| Assembled genome - Fasta file with chromosomes/scaffolds/contigs sequences     | X    |      |       |
| Gene annotation -  GTF file with gene annotations (UTRs,exons,CDS)     | X     |       |       |
| TE annotation -  GTF file with TE insertions     | X     |       |       |
| Reference transcripts - Fasta file with reference transcripts     |      |   X    |    X   |
| Reference TEs - Fasta with ref. TE insertions    |      |    X   |       |
| Dfam taxonomy OR fasta with ref. TE consensuses    |      |       |    X   |

---

## ChimeraTE genome-guided - Mode1 <a name="mode1"></a>
In the Mode 1, chimeric transcripts will be detected considering the genomic location of TE insertions and exons. Chimeras from this Mode can be classified as TE-initiated TE-exonized, and TE-terminated transcripts. Mode 1 does not detect chimeric transcripts derived from TE insertions absent from the reference genome that is provided. 
````
cd ChimeraTE/
python3 chimTE_mode1.py --help
````
````
ChimeraTE Mode 1: The genome-guided approach to detect chimeric transcripts with RNA-seq data.

Required arguments:
  --genome      Genome in fasta
  --input       Paired-end files and their respective group/replicate
  --project     Directory name with output data
  --te          GTF file containing TE information
  --gene        GTF file containing gene information
  --strand      Define the strandness direction of the RNA-seq. Two options:
                "rf-stranded" OR "fwd-stranded"

Optional arguments:
  --window      Upstream and downstream window size (default = 3000)
  --replicate   Minimum recurrency of chimeric transcripts between RNA-seq
                replicates (default 2)
  --coverage    Minimum coverage (mean between replicates default 2 for
                chimeric transcripts detection)
  --fpkm        Minimum fpkm to consider a gene as expressed (default 1)
  --threads     Number of threads (default 6)
  --overlap     Minimum overlap between chimeric reads and TE insertions (default 0.50)
  --index       Absolute path to pre-existing STAR index
````
---

### Prepare your data for Mode 1! <a name="prep_data"></a>
#### Input table
The input tab-delimited table provided with ````--input```` must have a specific format:
First column: Mate 1 from the paired-end data
Second column: Mate 2 from the paired-end data
Third column: Replicate/group name

| mate1 | mate2 | rep | 
| -------- | -------- | -------- |
| /home/user/ChimeraTE/mate1_control1.fastq.gz | /home/user/ChimeraTE/mate2_control1.fastq.gz | rep1 |
| /home/user/ChimeraTE/mate1_control2.fastq.gz | /home/user/ChimeraTE/mate2_control2.fastq.gz | rep2 |
| /home/user/ChimeraTE/mate1_control3.fastq.gz | /home/user/ChimeraTE/mate2_control3.fastq.gz | rep3 |

The header must be absent, as it follows in the example ````--input```` table at example_data/mode1/input_example.tsv

#### GTF for TEs
Usually, the coordinates for TE insertions is given as the .out file from RepeatMasker in many databases. If you already have a .out file from RepeatMasker, you can convert it to .gtf on Linux with:
````
tail -n +4 RMfile.out | egrep -v 'Satellite|Simple_repeat|rRNA|Low_complexity|RNA|ARTEFACT' | awk -v OFS='\t' '{Sense=$9;sub(/C/,"-",Sense);$9=Sense;print $5,"RepeatMasker","similarity",$6,$7,$2,$9,".",$10}' > RMfile.gtf
````
If you don't have the .out file for your genome assembly, check it out the util section.

---
 
### Example Data Mode 1 <a name="example_m1"></a>
After installation, you can run ChimeraTE with the example data from the sampled RNA-seq from *D. melanogaster* used in our paper.

````
#Do not forget to activate your conda environment:
conda activate chimeraTE
````
````
#One-line
python3 chimTE_mode1.py --genome example_data/mode1/dmel_genome_sample.fa --input example_data/mode1/input_mode1.tsv --project example_mode1 --te example_data/mode1/dmel_TEs_sample.gtf --gene example_data/mode1/dmel_genes_sample.gtf --strand rf-stranded

#Multi-line
python3 chimTE_mode1.py --genome example_data/mode1/dmel_genome_sample.fa \
--input example_data/mode1/input_mode1.tsv \
--project example_mode1 \
--te example_data/mode1/dmel_5TEs_sample.gtf \
--gene example_data/mode1/dmel_5genes_sample.gtf \
--strand rf-stranded
````

If you have more than 6 threads available on your machine, you can use ```--threads``` to speed up the process. 

---

### Output Mode 1 <a name="output_m1"></a>
The output files can be found at ```ChimeraTE/projects/$your_project_name```. For instance, for the example data, you can find the output at ```ChimeraTE/projects/example_mode1```. Inside this directory, you might found 3 tables:

 -  TE-initiated_final.ct
 -  TE-exonized_final.ct
 -  TE-terminated_final.ct

These tables contain the chimeric transcripts list with the location of genes and TE insertions generating chimeras, as well as their corresponding coverage of chimeric reads (support). At the 7th column of ```TE-exonized_final.ct```, you can find the position of the TE within the gene region (Embedded, Intronic, or Overlapped). As it follows in the example below:

=========================> TE-initiated_final.ct <=========================

| gene_id | gene_strand | gene_pos | TE_id | TE_strand | TE_pos | chim_reads | 
| -------- | -------- | -------- | -------- | -------- | -------- |  -------- |
| FBgn0031188 | - | X_RaGOO:21340686-21343686 | S2 | + | X_RaGOO:21341507-21342141 | 11.5 | 


=========================> TE-exonized_final.ct <=========================

| gene_id | gene_strand | gene_pos | TE_id | TE_strand | TE_pos | exonized_type | chim_reads | 
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| FBgn0285926 | - | X_RaGOO:10476773-10513188 | roo | - | X_RaGOO:10485868-10485985 | Embedded | 63.5 |
| FBgn0052000 | + | 4_RaGOO:126456-137357 | 1360 | + | 4_RaGOO:133965-134061 | Overlapped | 4.5 |
| FBgn0039923 | - | 4_RaGOO:761931-772400 | FB | - | 4_RaGOO:769101-769563 | Intronic | 91.0 |


=========================> TE-terminated_final.ct <=========================

| gene_id | gene_strand | gene_pos | TE_id | TE_strand | TE_pos | chim_reads | 
| -------- | -------- | -------- | -------- | -------- | -------- |  -------- |
| FBgn0011747 | - | 4_RaGOO:106334-1093346 | G5 | - | 4_RaGOO:109144-109334 | 5.0 |


## ChimeraTE genome-blinded - Mode 2 <a name="mode2"></a>

Mode 2 is designed to identify chimeric transcripts without the reference genome, with the prediction of chimeras from fixed and polymorphic TEs. In Mode 2, two alignments with stranded RNA-seq reads are performed: (1) against transcripts; (2) against TE insertions. From these alignments, all reads supporting chimeric transcripts (chimeric reads) will be computed. These reads are thise ones that have different singleton mates from the same read pairs splitted between transcripts and TEs, or those that have concordant alignment in one of the alignments, but singleton aligned reads in the other. There is also an option to perform de novo transcriptome assembly with ```--assembly``` parameter. Such additional analysis will analyze whether gene transcripts contain TE-derived sequences.

```
cd ChimeraTE/
python3 chimTE_mode2.py --help
```
```
ChimeraTE Mode 2: The genome-blinded approach to detect chimeric transcripts with RNA-seq data.

Required arguments:
  --input         Paired-end files and their respective group/replicate
  --project       Directory name with output data
  --te            Fasta file containing TE information
  --transcripts   Fasta file containing gene information
  --strand        Define the strandness direction of the RNA-seq. Two options:
                  "rf-stranded" OR "fwd-stranded"

Optional arguments:
  --coverage      Minimum coverage (mean between replicates default 2 for
                  chimeric transcripts detection)
  --fpkm          Minimum fpkm to consider a gene as expressed (default = 1)
  --replicate     Minimum recurrency of chimeric transcripts between RNA-seq
                  replicates (default = 2)
  --threads       Number of threads (default = 6)
  --assembly      Search for chimeric transcript with transcriptome assembly
                  with Trinity
  --ref_TEs       "species" database used by RepeatMasker (flies, human,
                  mouse, arabidopsis; or a built TE library in fasta format)
  --ram           Ram memory in Gbytes 
                  (default = 8)
  --overlap       Minimum overlap between chimeric reads and TE insertions
                  (default 0.50)
  --TE_length     Minimum TE length to keep it from RepeatMasker output
                  (default = 80bp)
  --identity      Minimum identity between de novo assembled transcripts and
                  reference transcripts (default = 80)
```

## Prepare your data for Mode 2 <a name="prep_mode2"></a>
Despite the format of the input files are simple fastas, altogether with paired-end RNA-seq reads, the sequence IDs for transcripts and TEs must be in a specific pattern. In order make it easier to generate these formats, we provide ```util``` scripts to manage your data.

#### 1. Reference transcripts (.fasta)
  - In order to run ChimeraTE correctly, this fasta file **must** have a specific header pattern. All IDs have be composed firstly by the isoform ID, followed by the gene name.  For instance, in _D. melanogaster_, the gene FBgn0263977 has two transcripts:<br />Tim17b-RA_FBgn0263977<br />Tim17b-RB_FBgn0263977
  - Note that headers "Tim17b-RA" and "Tim17b-RB" have isoform ID separated from gene name by "_".  This is not a usual ID format, thefore we have developed auxiliary scripts ($FOLDER/ChimeraTE/util/) to convert native ID formats to ChimeraTE format. 
    - *transcripts_IDs_NCBI.sh*     (native IDs from NCBI to the ChimeraTE format)
    - *transcripts_IDs_ensembl.sh*  (native IDs from ENSEMBL to the ChimeraTE format)
    - *transcripts_IDs_FLYBASE.sh*  (native IDs from FLYBASE to the ChimeraTE format)


### Example Data Mode 2 <a name="example_m2"></a>
After installation, you can run ChimeraTE Mode 2 with the example data from the sampled RNA-seq from *D. melanogaster* used in our paper.

````
#Do not forget to activate your conda environment:
conda activate chimeraTE
````

````
#One-line
python3 chimTE_mode2.py --input example_data/mode2/input_mode2.tsv --project example_mode2 --te example_data/mode2/dmel-sampled_TE-copies.fa --transcripts example_data/mode2/dmel-sampled_transcripts.fa --strand rf-stranded --assembly

#Multi-line
python3 chimTE_mode2.py --input example_data/mode2/input_mode2.tsv \
--project example_mode2 \
--te example_data/mode2/dmel-sampled_TE-copies.fa\
 --transcripts example_data/mode2/dmel-sampled_transcripts.fa \
 --strand rf-stranded \
--assembly
````
Mode 2 will run with 8 threads and 8Gb of RAM memory, but you can speed up the analysis by increasing this values with ```--threads``` and ```--ram```, respectively.

**NOTE**: If you are not working with *Drosophila* data, do not forget to change ```--ref_TEs```parameter, providing a Dfam taxonomy level to use with RepeatMasker, or a fasta with TE consensuses.


### Output Mode 2 <a name="output_m2"></a>
The output files can be found at ```ChimeraTE/projects/$your_project_name```. For instance, for the example data, you can find the output at ```ChimeraTE/projects/example_mode2```. Inside this directory, you might found 3 tables:

 -  chimreads_evidence_FINAL.tsv
<br />In the "chimreads_evidence" table, you will find chimeric transcripts supported **only** by paired-end reads that have mapped in both transcripts and TE sequences (singletons and concordant/singleton - Check manuscripts's methods).
 -  transcriptome_evidence_FINAL.tsv
<br />In the "transcriptome_evidence" table, you will find chimeras supported **only** by the transcripme assembly method (if you have activated ```--assembly```option). This table will provide you the gene, TE family, and the respective assembled transcript ID for which a TE sequence was found.
 -  double_evidence_FINAL.tsv
<br />Finally, "double_evidence" is the list of chimeras for which both previous methods have predicted the same chimera (strong evidence!), containing all information from both previous tables.


=========================> chimreads_evidence_FINAL.tsv <=========================

| gene_id | TE_family | chim_reads | transcript_ID | transcript_FPKM | 
| -------- | -------- | -------- | -------- | -------- |
| FBgn0058160 | DNAREP1 | 60.0 | CG40160-RH_FBgn0058160 | 62177.475 |

=========================> transcriptome_evidence_FINAL.tsv <=========================

| gene_id | TE_family | transcript_ID | Trinity_transcripts | Identity_transcripts | trinity_length | ref_transcript_length | match_length | chim_reads | 
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| FBgn0286778 | HMSBEAGLE_I |	CG46385-RA | TRINITY_DN87_c0_g1_i1; TRINITY_DN88_c0_g1_i1 |	97.992 |	741.5 | 5129.0 | 732.0 | 31.0 |


=========================> double_evidence_FINAL.tsv <=========================

| gene_id | TE_family | chim_reads | masked_family | chim_reads_masked | ref_transcript_FPKM | Trinity_transcripts | Identity_transcripts | trinity_length | ref_transcript_length | match_length | ref_transcript_IDs | 
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| FBgn0001169 | ROO | 32.0 | ROO_I | 4.0 | 3011.5781 | TRINITY_DN13_c0_g1_i3; TRINITY_DN13_c0_g1_i2 |	100.0 | 604.5 | 4069.5 | 603.5 |	H-RD; H-RB	H-RD_FBgn0001169; H-RB_FBgn0001169; H-RA_FBgn0001169 |

---

## Please cite us!
Daniel S Oliveira, Marie Fablet, Anaïs Larue, Agnès Vallier, Claudia M A Carareto, Rita Rebollo, Cristina Vieira. ChimeraTE: A pipeline to detect chimeric transcripts derived from genes and transposable elements. Nucleic and Acids Research, 2023. https://doi.org/10.1093/nar/gkad671

## Development and help
To report bugs and give us suggestions, you can open an [issue](https://github.com/OliveiraDS-hub/ChimeraTE/issues) on the github repository.

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/80x15.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
