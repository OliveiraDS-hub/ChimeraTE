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
  
### Util scripts

#### Masking the genome with RepeatMasker (output: proper .gtf to use as mode1 input)
````
bash util/masking.sh [--genome <genome.fa>] [--ref_TEs <flies/mouse/human>] [--out <output_file>] [options]

#Mandatory arguments:

  --genome    file with genome (.fasta)

  --ref_TEs   "species" database used by RepeatMasker (flies, human, mouse, arabidopsis...)

  --out       output file

#Optional arguments:

  --threads   Number of threads, (default: 6)

  --dist      Distance in nt between LTR and Internal regions, as well as fragments from the
              same TE family that will be merged, (default: 50)
````
  
#### .out file from RepeatMasker to fasta (output: proper fasta to use as mode2 input)
````
bash util/rmout2fasta.sh [--genome <genome.fa>] [--rm <repeatmasker.out>] [--out <output_file>
  #Mandatory arguments:

  --genome    file with genome (.fasta)
  --rm   file from RepeatMasker (.out)
  --out   output with TE insertions (.fasta)
````
  
### ChimeraTE genome-guided (mode1)
````
bash ChimeraTE_mode1.sh --help
````

````
USAGE:

-One-replicate:

ChimTE-mode1.sh [--mate1 <mate1.fastq.gz>] [--mate2 <mate2.fastq.gz>] [--genome <genome.fasta>] [--te <TE_insertions.gtf or TE_insertions.bed>] [--gene <gene_annotation.gtf>] [--project <project_name>]

#Mandatory arguments:

   --mate1			mate 1 from paired-end reads

   --mate2			mate 2 from paired-end reads

   --genome			genome sequence .fa

   --te				GTF file with TE coordinates

   --gene			GTF file with genes coordinates

   --project			project name

#Optional arguments:

   --window			Upstream and downstream window size (default = 3000)

   --overlap  			Minimum overlap between chimeric reads and TE insertions

   --utr 			It must be used if your gene annotation (-a | --gene) has UTR regions (default = off)

   --fpkm   			Minimum fpkm to consider a gene as expressed (default = 1)

   --threads			Number of threads (default:6)
 ````
 
Explain the pipeline here

### ChimeraTE de novo approach (mode2)

This mode is going to perform two alignment with stranded RNA-seq reads: (1) against transcripts; (2) against TE insertions. From these alignments, all reads supporting chimeric transcripts (chimeric reads) will be computed. These reads are thise ones that have different singleton mates from the same read pairs splitted between transcripts and TEs, or those that have concordant alignment in one of the alignments, but singleton aligned reads in the other.


In order to run this mode, despite the format of the input files are simple fastas, the sequence IDs must be in a specific pattern. 
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
````
bash ChimeraTE_mode2.sh --help
````

````
USAGE:

-One-replicate:

ChimTE-mode2.sh [--mate1 <mate1.fastq.gz>] [--mate2 <mate2.fastq.gz>] [--te <TE_insertions.fa>] [--transcripts <transcripts.fa>] [--stranded <rf-stranded or fr-stranded>] [--project <project_name>] [options]

-Multi-replicates:

ChimTE-mode2.sh [--mate1 <mate1_replicate1.fastq.gz,mate1_replicate2.fastq.gz>] [--mate2 <mate2_replicate1.fastq.gz,mate2_replicate2.fastq.gz>] [--te <TE_insertions.fa>] [--transcripts <transcripts.fa>] [--stranded <rf-stranded or fr-stranded>] [--project <project_name>] [options]


#Mandatory arguments:

  --mate1 			mate 1 from paired-end reads

  --mate2 			mate 2 from paired-end reads

  --te				TE insertions (.fasta)

  --transcripts			transcripts (.fasta)

  --stranded			Select "rf-stranded" if your reads are reverse->forward; or "fr-stranded" if they are forward->reverse

  --project			project name

#Optional arguments:

  --cutoff			Minimum chimeric reads as support (default: 2)

  --threads			Number of threads, (default: 6)
````

### ChimeraTE transcriptome-guided (mode3)

````
bash ChimeraTE_mode3.sh --help
````

````
USAGE:

-One-replicate:

ChimTE-transcripts.sh [--mate1 <mate1.fastq.gz>] [--mate2 <mate2.fastq.gz>] [--stranded <RF or FR>] [--transcripts <transcripts.fa>] [--ref_TEs <taxonomy or library_TEs.fa>] [options]

-Multi-replicates:

ChimTE-transcripts.sh [--mate1 <mate1_replicate1.fastq.gz,mate1_replicate2.fastq.gz>] [--mate2 <mate2_replicate1.fastq.gz,mate2_replicate2.fastq.gz>] [--stranded <RF or FR>] [--transcripts <transcripts.fa>] [--ref_TEs <taxonomy or library_TEs.fa>] [options]


#Mandatory arguments:

    --mate1			mate 1 from paired-end reads

    --mate2			mate 2 from paired-end reads

    --stranded		  	indicates the order of the stranded RNA-seq (RF - reverse/forward; FR - forward/reverse)

    --transcripts		fasta file with reference transcripts

    --ref_TEs			"species" database used by RepeatMasker (flies, human, mouse, arabidopsis);
				or a built TE library in fasta format

#Optional arguments:

    --output		  	Output directory

    --threads               	Number of threads, (default: 6)

    --ram			RAM memory (default: 8 Gbytes)

    --TE_length                 Minimum TE length to keep it from RepeatMasker output (default: 80bp)

    --min_length          	Minimum length with homology between de novo assembled transcripts and reference transcripts (default: 80%)

    --overlap              	Minmum overlap (0.1 to 1) between read length and TE insertion (default: 0.5)

    --min_TPM             	Minimum TPM expression (default: 1)

    -h, --help                  Show this help message and exit
````

## License
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/80x15.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
