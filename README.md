# The pipeline

ChimeraTE is a pipeline to detect chimeric transcripts derived from genes and transposable elements (TEs). It has two usage Modes:

Mode 1 chimeric transcripts detection based upon exons and TE copies positions in the genome sequence; 

Mode 2 chimeric transcripts detection regardless the genomic position, allowing the detection of chimeras from TEs that are not present in the referece genome, but with less sensivity.

It has been tested in Linux machines, Ubuntu 18.04 and 20.04.

1. [Install](#installation)
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

---

## Install <a name="installation"></a>
The installation may be easily done with conda. If you don't have conda installed in your machine, please follow [this tutorial](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html).

Once you have installed conda, all dependencies to run ChimeraTE can be easily installed in a new conda environment by using the chimeraTE.yml file:
````
#Download repository from github
git clone https://github.com/oliveirads-bioinfo/ChimeraTE.git

#Go to ChimeraTE's folder
cd $FOLDER/ChimeraTE

#create chimeraTE environment with all dependencies
conda env create -f chimeraTE.yml

#activate the new environment
conda activate chimeraTE

#Give write permissions
chmod a+x $FOLDER/ChimeraTE/scripts/*

````
---

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
| TE consensuses - Fasta file with reference TE consensuses    |      |       |    X   |

---

## ChimeraTE genome-guided (mode1) <a name="mode1"></a>
````
cd $FOLDER/ChimeraTE/
bash ChimeraTE-mode1.sh --help
````

### USAGE <a name="usage"></a>
````
ChimTE-mode1.sh   [--mate1 <replicate1_R1.fastq.gz,replicate2_R1.fastq.gz,replicate3_R1.fastq>]
                  [--mate2 <replicate1_R2.fastq.gz,replicate2_R2.fastq.gz,replicate3_R2.fastq>]
                  [--genome <genome.fasta>]
                  [--te <TE_insertions.gtf>]
                  [--gene <gene_annotation.gtf>]
                  [--project <project_name>]

#Mandatory arguments:

   --mate1                 FASTQ paired-end R1 (replicates separated by ",")

   --mate2                 FASTQ paired-end R2 (replicates separated by ",")

   --genome                FASTA genome sequence

   --te                    GTF/GFF with TE coordinates

   --gene                  GTF/GFF with genes coordinates

   --project               project name (it's the name of the directory that will be created inside projects/)

#Optional arguments:

   --window                Upstream and downstream window size (default = 3000)

   --overlap               Minimum overlap between chimeric reads and TE insertions (default = 0.5 -50%-)

   --utr                   It must be used if your gene annotation (-a | --gene) has UTR regions (default = off)

   --fpkm                  Minimum fpkm to consider a gene as expressed (default = 1)

   --coverage              Minimum coverage for chimeric transcripts detection (default = 2)

   --replicate	           Minimum recurrency of chimeric transcripts between RNA-seq replicates (default = 2)

   --threads               Number of threads (default = 6)
````
---

### Prepare your data for Mode 1! <a name="prep_data"></a>

Stranded paired-end libraries, genome fasta file and gene annotation are common files that can be found in databases for many species. However, TE annotation may not be found publicly available to non-reference species/strains. 

Then, if you don't have TE annotation for your assembled genome, you can use our *masking.sh* util script to mask it with RepeatMasker and One Code to Find Them All.
The genome is masked with TE consensus available in Dfam 3.6, by providing the lineage with ````--ref_TEs```` or a fasta file with a TE consensus library.

````
bash util/masking.sh [--genome <genome.fa>] [--ref_TEs <flies/mouse/human OR TE_library.fa>] [--out <output_file>] [options]
````
| Parameter | Description |
| -------- | -------- |
| --genome     |  Fasta file with chromosomes/scaffolds/contigs sequences |
| --ref_TEs     |  Fasta file with TE consensuses <br />OR<br /> RepeatMasker Dfam reference (i.e.:flies,arabidopsis,human)|
| --out     |  output file (gtf format) |
| --threads     |  Number of threads to run RepeatMasker (default: 6) |
| --dist     | Distance in nt between LTR and Internal regions, as well as <br />fragments from the same TE family that will be merged, (default: 50)  |

This script will provide you a gtf file with TE annotation to run ChimeraTE Mode 1.

If you already have a .out file from RepeatMasker, you can convert it to .gtf with:
````
tail -n +4 RMfile.out | awk -v OFS='\t' '{Sense=$9;sub(/C/,"-",Sense);$9=Sense;print $5,"RepeatMasker","similarity",$6,$7,$2,$9,".",$10}' > RMfile.gtf
````

---
 
### Example Data Mode 1 <a name="example_m1"></a>
After installation, you can run ChimeraTE with the example data from the sampled RNA-seq from *D. melanogaster* used in our paper.

````
cd $FOLDER/ChimeraTE/example_data/mode1
gunzip *
cd ../../

bash ChimeraTE-mode1.sh --mate1 example_data/mode1/sample1_R1.fq,example_data/mode1/sample2_R1.fq \
                        --mate2 example_data/mode1/sample1_R2.fq,example_data/mode1/sample2_R2.fq \
                        --genome example_data/mode1/dmel-all-chromosome-r6.46.fasta \
                        --te example_data/mode1/dmel-all-chromosome-r6.46_RM_final.gtf \
                        --gene example_data/mode1/dmel-all-r6.46.gtf \
                        --project sampling-mode1 \
                        --utr
````
---

### Output Mode 1 <a name="out_m1"></a>
blabla

---

## ChimeraTE genome-blinded (mode2) <a name="mode2"></a>

This mode is going to perform two alignment with stranded RNA-seq reads: (1) against transcripts; (2) against TE insertions. From these alignments, all reads supporting chimeric transcripts (chimeric reads) will be computed. These reads are thise ones that have different singleton mates from the same read pairs splitted between transcripts and TEs, or those that have concordant alignment in one of the alignments, but singleton aligned reads in the other.

  
````
cd $FOLDER/ChimeraTE/
bash ChimeraTE_mode2.sh --help
````

### USAGE <a name="usage_m2"></a>
````
ChimTE-mode2.sh [--mate1 <replicate1_R1.fastq.gz,replicate2_R1.fastq.gz,replicate3_R1.fastq>]
                [--mate2 <replicate1_R2.fastq.gz,replicate2_R2.fastq.gz,replicate3_R2.fastq>] 
                [--te <TE_insertions.fa>] 
                [--transcripts <transcripts.fa>] 
                [--stranded <rf-stranded or fr-stranded>] 
                [--project <project_name>]

#Mandatory arguments:

  --mate1                 FASTQ paired-end R1 (replicates separated by ",")

  --mate2                 FASTQ paired-end R2 (replicates separated by ",")

  --te                    TE insertions (fasta)

  --transcripts           transcripts (fasta)

  --stranded              Select "rf-stranded" if your reads are reverse->forward; or "fr-stranded" if they are forward->reverse

  --project               project name

#Optional arguments:

  --fpkm                  Minimum FPKM expression (default = 1)

  --coverage              Minimum chimeric reads as support (default = 2)

  --replicate	          Minimum recurrency of chimeric transcripts between RNA-seq replicates (default = 2)

  --threads               Number of threads, (default = 6)

  --assembly              Perform transcripts assembly -HIGH time consuming- (default = deactivated)
       |
       |
       -------> Mandatory if --assembly is activated:
       |        --ref_TEs             "species" database used by RepeatMasker (flies, human, mouse, arabidopsis);
       |                              or a built TE library in fasta format
       |
       |
       -------> Optional if --assembly is activated:
                --ram                 RAM memory (default: 8 Gbytes)
                --overlap             Minmum overlap (0.1 to 1) between read length and TE insertion (default = 0.5)
                --TE_length           Minimum TE length to keep it from RepeatMasker output (default = 80bp)
                --min_length          Minimum identity between de novo assembled transcripts and reference transcripts (default = 80%)
````

---

### Prepare your data to Mode 2 <a name="prep_data_m2"></a>

Despite the format of the input files are simple fastas, the sequence IDs must be in a specific pattern. 
  
#### 1. Reference transcripts (.fasta)
  - In order to run ChimeraTE correctly, this fasta file **must** have a specific header pattern. All IDs must be composed firstly by the isoform ID, followed by the gene name. For instance, in _D. melanogaster_, the gene FBgn0263977 has two transcripts:<br />Tim17b-RA_FBgn0263977<br />Tim17b-RB_FBgn0263977
  - Note that headers "Tim17b-RA" and "Tim17b-RB" have isoform ID separated from gene name by "_".  This is not a usual ID format, thefore we have developed auxiliary scripts ($FOLDER/ChimeraTE/util/) to convert native ID formats to ChimeraTE format. 
    - *transcripts_IDs_NCBI.sh*     (native IDs from NCBI to the ChimeraTE format)
    - *transcripts_IDs_ensembl.sh*  (native IDs from ENSEMBL to the ChimeraTE format)
    - *transcripts_IDs_FLYBASE.sh*  (native IDs from FLYBASE to the ChimeraTE format)

````
cd $FOLDER/ChimeraTE/util
bash transcripts_IDs_[NCBI-ENSEMBL-FLYBASE].sh --help
````
````
Conversion of [NCBI-ENSEMBL-FLYBASE] native transcript IDs to the ChimeraTE Mode 2 format

#Mandatory arguments:

  --transcripts     transcripts downloaded from [NCBI-ENSEMBL-FLYBASE] (.fasta)

  --out   output file name (.fasta)
````
We provide here the corrected IDs for *D. melanogaster*, human (hg38), mouse and *A. thaliana*. 

#### 2. Reference TE insertions (.fasta)
- The TE (.fasta) file used by Mode 2 must have only TE insertions. Be sure that they do not contains any Satellites or Tandem repeats.
In addition, the reference TE insertions **must** have only the TE family in the headers. For instance, if *D. melanogaster* genome has ~4.000 DNAREP-1 TE insertions, all of them must have the header as ">DNAREP-1". If you have the .out file from RepeatMasker, you can generate the fasta file with the proper headers to run ChimeraTE Mode 2 with *rmout2fasta.sh* util script.
    - *rmout2fasta.sh*      (RepeatMasker .out file to .fasta file with ChimeraTE Mode 2 IDs)

````
cd $FOLDER/ChimeraTE/util
bash util/rmout2fasta.sh --help 
````
````
rmout2fasta.sh [--genome <genome.fa>] [--rm <repeatmasker.out>] [--out <output_file>]
  #Mandatory arguments:

  --genome    file with genome (.fasta)
  --rm   file from RepeatMasker (.out)
  --out   output with TE insertions (.fasta)
````
We provide here the corrected fasta file with all headers formatted for _D. melanogaster_, human (hg38), mouse (mmX) and _A. thaliana_. 

---

### Example data Mode 2 <a name="example_m2"></a>
````
cd $FOLDER/ChimeraTE/example_data/mode1
gunzip *
cd ../../

bash ChimeraTE-mode2.sh --mate1 example_data/mode1/sample1_R1.fq,example_data/mode1/sample2_R1.fq \
                        --mate2 example_data/mode1/sample1_R1.fq,example_data/mode1/sample2_R1.fq \
                        --te data/data_sampling-MODE2/dmel-all-chromosome-r6.43_RM_final.fasta \
                        --transcripts data/data_sampling-MODE2/dmel-all-transcripts-clean.fa \
                        --stranded rf-stranded \
                        --project sampling-mode2 \
                        --assembly \
                        --ref_TEs util/bergman_cons.fa \
````
---

### Output Mode 2 <a name="out_m2"></a>
blabla

---

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/80x15.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
