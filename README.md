# The pipeline

ChimeraTE is a pipeline to detect chimeric transcripts derived from genes and transposable elements (TEs). It has two running Modes:

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
git clone https://github.com/OliveiraDS-hub/ChimeraTE.git

#Go to ChimeraTE's folder
cd $FOLDER/ChimeraTE

#create chimeraTE environment with all dependencies
conda env create -f chimeraTE.yml

#activate the new environment
conda activate chimeraTE

#Give write permissions
chmod +x scripts/mode1/* scripts/mode2/*
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
In the Mode 1, chimeric transcripts will be detected considering the genomic location of TE insertions and exons. Chimeras from this Mode can be classified as TE-initiated upstream, TE-initiated 5’UTR, TE-exonized, TE-terminated 3’UTR and TE-terminated downstream. In addition, results from Mode 1 can be visualized in genome browsers, which allows a manual curation of chimeric transcripts in the reference genome. Mode 1 does not detect chimeric transcripts derived from TE insertions absent from the reference genome that is provided. 
````
cd $FOLDER/ChimeraTE/
bash ChimeraTE-mode1.sh --help
````
### USAGE <a name="usage"></a>
````
ChimTE-mode1.sh   [--mate1 <replicate1_R1.fastq.gz,replicate2_R1.fastq.gz,replicate3_R1.fastq>]
                  [--mate2 <replicate1_R2.fastq.gz,replicate2_R2.fastq.gz,replicate3_R2.fastq>]
                  [--strand <rf-stranded> OR <fr-stranded>]
                  [--genome <genome.fasta>]
                  [--te <TE_insertions.gtf>]
                  [--gene <gene_annotation.gtf>]
                  [--project <project_name>]

#Mandatory arguments:

   --mate1                 FASTQ paired-end R1

   --mate2                 FASTQ paired-end R2

   --genome                FASTA genome sequence

   --te                    GTF/GFF with TE coordinates

   --gene                  GTF/GFF with genes coordinates

   --strand                Select "rf-stranded" if your reads are reverse->forward; or "fr-stranded" if they are forward->reverse

   --project               project name -spaces are not allowed- (it's the name of the directory that will be created inside projects/)

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
   - If you don't have TE annotation for your assembled genome, you can use the util script *masking.sh*.The genome is masked with TE consensus available in Dfam 3.6, by providing the lineage with ````--ref_TEs```` or a fasta file with a TE consensus library. This script will provide you a gtf file with TE annotation to run ChimeraTE Mode 1.

````
cd $FOLDER/ChimeraTE/util
bash masking.sh --help
````
````
#Mandatory arguments:

  --genome    file with genome (.fasta)

  --ref_TEs   species database used by RepeatMasker ("flies", "human", "mouse", "arabidopsis"...);
  	      or a file with your custom TE library (.fasta)

#Optional arguments:

  --out       output directory (if not provided, the output will be created in the current folder)

  --threads   Number of threads, (default = 6)

  --dist      Distance in nt between LTR and Internal regions, as well as fragments from the
              same TE family that will be merged by One Code to Find them all, (default = 50)
````

   - If you already have a .out file from RepeatMasker, you can convert it to .gtf with:
````
tail -n +4 RMfile.out | awk -v OFS='\t' '{Sense=$9;sub(/C/,"-",Sense);$9=Sense;print $5,"RepeatMasker","similarity",$6,$7,$2,$9,".",$10}' | egrep -v 'Satellite|Simple_repeat|rRNA|Low_complexity|RNA|ARTEFACT' > RMfile.gtf
````

---
 
### Example Data Mode 1 <a name="example_m1"></a>
After installation, you can run ChimeraTE with the example data from the sampled RNA-seq from *D. melanogaster* used in our paper.

````
cd $FOLDER/ChimeraTE/example_data/mode1
gunzip *
cd ../../
                        
bash ChimeraTE-mode1.sh --mate1 example_data/data_sampling-MODE1/sample1_R1.fq.gz,example_data/data_sampling-MODE1/sample2_R1.fq.gz \
--mate2 example_data/data_sampling-MODE1/sample1_R2.fq.gz,example_data/data_sampling-MODE1/sample2_R2.fq.gz \
--genome example_data/data_sampling-MODE1/dmel-chrX.fa \
--te example_data/data_sampling-MODE1/dmel_sample-TEs-chrX.gtf \
--gene example_data/data_sampling-MODE1/dmel_sample-transcripts-chrX.gtf \
--strand rf-stranded \
--project example_data-Mode1 \
--utr
````
---

### Output Mode 1 <a name="out_m1"></a>
The output files can be found at ```$DIR/ChimeraTE/projects/$your_project_name```. For instance, for the example data, you can find the output at ```$DIR/ChimeraTE/projects/sampling-mode1```. Inside this directory, you might found 5 tables:
   
   - TE-initiated-UP_final.ct
   - TE-initiated-5UTR_final.ct
   - TE-exonized_final.ct
   - TE-terminated-3UTR_final.ct
   - TE-terminated-DOWN_final.ct

Each table have a list of chimeric transcripts based on chimeric reads evidence. For TEs located upstream/downstream and with chimeric reads with exons (UTRs/CDSs), they can be found in ```TE-initiated-UP_final.ct``` and ```TE-terminated-DOWN_final.ct```, respetively. For TEs located inside the gene region, they are separated into three categories: (1) Evidence of chimeric reads between TE and 5'UTR in ```TE-initiated-5UTR_final.ct```; (2) chimeric reads between TE and CDS region in ```TE-exonized_final.ct``` and (3) chimeric reads between TE and 3' UTR in ```TE-terminated-3UTR_final.ct```.
If in your dataset one of these tables is missing, it means that there were no chimeric transcripts detected for these category.

Checking the ```TE-exonized_final.ct``` from example data:

    head $DIR/ChimeraTE/projects/sampling-mode1/TE-exonized_final.ct

You should find at the first lines:
 
| gene_id | gene_strand | TE_family | TE_strand |  TE_position | Exon-position | chimeric_reads |
| -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| FBgn0001169 | +    | roo  | +   | 3R:20622645-20622806 | 3R:20622584-20623298 | 49.0 | 
| FBgn0003028 | +    | roo  | +   | X:5066743-5066982 | X:5066073-5067781 | 61.0 | 


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
````
To run it with the transcriptome assembly option:
````
bash ChimeraTE-mode2.sh --mate1 example_data/data_sampling-MODE2/sample1_R1.fq.gz,example_data/data_sampling-MODE2/sample2_R1.fq.gz \
--mate2 example_data/data_sampling-MODE2/sample1_R2.fq.gz,example_data/data_sampling-MODE2/sample2_R2.fq.gz \
--te example_data/data_sampling-MODE2/dmel-sampled_TE-copies.fa \
--transcripts example_data/data_sampling-MODE2/dmel-sampled_transcripts.fa \
--strand rf-stranded \
--project example_mode2 \
--threads 8 \
--assembly \
--ref_TEs example_data/data_sampling-MODE2/dmel-sampled_TEconsensus.fa
````
If you don't want to test the assembly option, just remove ````--assembly```` and ````--ref_TEs```` from the commandline.

---
### Output Mode 2 <a name="out_m2"></a>
The output files can be found at ```$DIR/ChimeraTE/projects/$your_project_name```. For instance, for the example data, you can find the output at ```$DIR/ChimeraTE/projects/sampling-mode2```. Inside this directory, you'll find 3 tables:
   - chimTE-final-chimreads.ct
   - chimTE-final-transcriptome.ct
   - chimTE-final-double-evidence.ct

````chimTE-final-chimreads.ct````:
This table contains all chimeric transcripts detected **only** through chimeric reads evidence. Chimeras listed on this table will be also found when ````--assembly```` option is disabled.

    head $DIR/ChimeraTE/projects/sampling-mode2/chimTE-final-chimreads.ct

You should find at the first lines:
| gene_id | TE_family | chimeric_reads |  ref_transcript | FPKM |
| -------- | -------- | -------- | -------- | -------- | 
| FBgn0000448 | ROO_I    | 6.0000  | Hr3-RA_FBgn0000448   | 79.3795 |
| FBgn0001169 | ROO_I    | 34.0000  | H-RA_FBgn0001169,H-RB_FBgn0001169,H-RD_FBgn0001169   | 367.3170 | 

````chimTE-final-transcriptome.ct````:
This table contains all chimeric transcripts detected **only** through the transcriptome assembly evidence.

    head $DIR/ChimeraTE/projects/sampling-mode2/chimTE-final-transcriptome.ct
 
You should find at the first lines:

| Trinity_isoform | Ref_transcript | Gene |  Identity | Chimeric_transcript_length | Ref_transcript_length | Match_length | TE_family |  Chimeric_reads |
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |
| TRINITY_DN83_c0_g1_i1 | Abl-RB    | FBgn0000017  | 99.800   | 5507 | 12564 | 6345 | 297 | 2.0000 |
| TRINITY_DN19_c0_g1_i2 | Agpat1-RB | FBgn0030421  | 98.718   | 477 | 3795 | 468 | 1360 | 3.0000 |

````chimTE-final-double-evidence.ct````:
This table contains all chimeric transcripts detected through **both** chimeric reads and transcriptome assembly evidences.

    head $DIR/ChimeraTE/projects/sampling-mode2/chimTE-final-double-evidence.ct
 
You should find at the first lines:
| Gene | TE_family | Chimeric_reads |  Transcripts | FPKM | Trinity_isoform | Identity | Chimeric_transcript_length |  Ref_Length | Match_length | Masked_TE_family | chim_reads_trinity|
| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- |  -------- | -------- | -------- |
| FBgn0004882 | ROO_I | 49.0000  | orb-RA_FBgn0004882,orb-RB_FBgn0004882,orb-RD_FBgn0004882,orb-RF_FBgn0004882,orb-RG_FBgn0004882   | 405.4890 | TRINITY_DN5_c1_g1_i1 | 100.000 | 474 | 4862 | 474 | roo | 63.0000 |
| FBgn0010215 | DNAREP1_DM | 11.0000  | alpha-Cat-RB_FBgn0010215 | 411.5380 | TRINITY_DN83_c0_g1_i1 | 99.582 | 479 | 3487 | 479 | INE-1 | 3.0000 |

---

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/80x15.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
