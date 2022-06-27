# ChimeraTE
## Install
The installation may be easily done with conda, with chimeraTE.yml file:
````
#create chimeraTE environment with all dependencies
conda env create -f chimeraTE.yml

#activate the new environment
conda activate chimeraTE
````
### Dependencies
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [seqtk](https://github.com/lh3/seqtk)
- [samtools](http://www.htslib.org/download/)
- [bedtools](https://github.com/arq5x/bedtools2/releases)
- [express](https://pachterlab.github.io/eXpress/overview.html#)
- [RSEM](https://github.com/deweylab/RSEM)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)

### Inputs:

  #### 1. Stranded paired-end RNA-seq
  - It's strongly recommended to use ChimeraTE with RNA-seq replicates.
  
  #### 2. Reference transcripts (.fasta)
  - In order to run ChimeraTE correctly, this fasta file **must** have a specific header pattern. All IDs must be composed firstly by the isoform ID, followed by the gene name. For instance, in _D. melanogaster_, the gene FBgn0263977 has two transcripts:
  Tim17b-RA_FBgn0263977
  Tim17b-RB_FBgn0263977
  
  Note that headers "Tim17b-RA" and "Tim17b-RB" have isoform ID separated from gene name by "_". 
  This is not a usual ID format, thefore we have developed auxiliary scripts (scripts/aux/) to convert native IDs the script _transcripts_IDs_NCBI.sh_ to transform the IDs from native NCBI format to the ChimeraTE format (see details in Manual). This script may be used if you are using a genome annotation from NCBI.

  In addition, we provide here the corrected IDs for _D. melanogaster_, human (hg38), mouse (mmX) and _A. thaliana_. 

  #### 3. Reference TE insertions (.fasta)

  - This .fasta file must have only TE insertions. Be sure that they do not contains any Satellites or Low complexity repeats.
  - The .fasta file with the reference TE insertions **must** have only the TE family in the headers. For instance, if _D. melanogaster_ genome has ~4.000 DNAREP-1 TE insertions, all of them must have the header as ">DNAREP-1".

  We provide here the corrected fasta file with all headers formatted for _D. melanogaster_, human (hg38), mouse (mmX) and _A. thaliana_. 
  
### Usage
````
bash chimeraTE_mdl1.sh --help
````

````
ChimeraTE - Usage commandline

   Required arguments:

        -1 | --mate1    	paired-end R1 (multiple replicates must be comma-separated)

        -2 | --mate2    	paired-end R2 (multiple replicates must be comma-separated)

        -g | --gene     	gene sequences .fa

        -t | --te       	TE sequences .fa

        -p | --project  	project name
  
        -s | --strandness	Select rf-stranded if your reads are reverse->forward; or fr-stranded if they are forward->reverse

   Optional arguments:

	      -c | --cutoff   	Minimum chimeric pairs
        
        # threads
        
        # min fpkm
 ````
 
 ### Example
 
 
