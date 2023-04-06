![image](https://github.com/OliveiraDS-hub/ChimeraTE/blob/main/image/logo.png)

**The previous bash version (ChimeraTE v1.0) is deprecated!!!**

A python version with several improvements (ChimeraTE v1.1) will be available soon.

---
# The pipeline
ChimeraTE is a pipeline to detect chimeric transcripts derived from genes and transposable elements (TEs). It has two running Modes:

- Mode 1 chimeric transcripts detection based upon exons and TE copies positions in the genome sequence; 

- Mode 2 chimeric transcripts detection regardless the genomic position, allowing the detection of chimeras from TEs that are not present in the referece genome, but with less sensivity.


1. [Install](#installation)
   1. [Conda](#conda)
   2. [Docker](#docker)
2. [Required data](#req_data)
3. [ChimeraTE Mode 1](#mode1)
    1. [Usage](#usage)
    2. [Preparing your data](#prep_data)
    3. [Example data](#example_m1)
    4. [Output](#out_m1)
4. [ChimeraTE Mode 2](#mode2)
    1. [Usage](#usage_m2)
    2. [Preparing your data](#prep_data_m2)
    3. [Example data](#example_m2)
    4. [Output](#out_m2)
5. [Util]

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

## Required data <a name="req_data"></a>
In order to run ChimeraTE, the following files are required according to the running Mode: 

| Data | Mode 1 | Mode 2 | Mode 2 --assembly |
| -------- | -------- | -------- | -------- |
| Stranded paired-end RNA-seq - Fastq files     | X    |   X   |    X   |
| Assembled genome - Fasta file with chromosomes/scaffolds/contigs sequences     | X    |      |       |
| Gene annotation -  GTF file with gene annotations (UTRs,exons,CDS)     | X     |       |       |
| TE annotation -  GTF file with TE insertions     | X     |       |       |
| Reference transcripts - Fasta file with reference transcripts     |      |   X    |    X   |
| Reference TEs - Fasta file with reference TE insertions    |      |   X    |    X   |

---

## ChimeraTE genome-guided (mode1) <a name="mode1"></a>
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
  --threads     Number of threads (default 6
  --overlap     Minimum overlap between chimeric reads and TE insertions
                (default 0.50)
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
tail -n +4 RMfile.out | awk -v OFS='\t' '{Sense=$9;sub(/C/,"-",Sense);$9=Sense;print $5,"RepeatMasker","similarity",$6,$7,$2,$9,".",$10}' | egrep -v 'Satellite|Simple_repeat|rRNA|Low_complexity|RNA|ARTEFACT' > RMfile.gtf
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
python3 chimTE_mode1.py --genome example_data/dmel-chrX.fa \
--input example_data/input_example.tsv \
--project example_test \
--te example_data/dmel_sample-TEs-chrX.gtf \
--gene example_data/dmel_sample-transcripts-chrX.gtf \
--strand rf-stranded
````
---

## Read more about ChimeraTE!
Check it out the pre-print paper on BioRchive!<br />ChimeraTE: A pipeline to detect chimeric transcripts derived from genes and transposable elements. https://doi.org/10.1101/2022.09.05.505575<br />

## Development and help
To report bugs and give us suggestions, you can open an [issue](https://github.com/OliveiraDS-hub/ChimeraTE/issues) on the github repository. Feel free also to contact me by e-mail (daniel.sdo2015@gmail.com) or twitter (https://twitter.com/Oliveira_D_S)

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/80x15.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
