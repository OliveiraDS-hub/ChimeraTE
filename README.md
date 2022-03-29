# ChimeraTE


## Install
### Dependencies
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [seqtk](https://github.com/lh3/seqtk)
- [samtools](http://www.htslib.org/download/)
- [bedtools](https://github.com/arq5x/bedtools2/releases)
- [express](https://pachterlab.github.io/eXpress/overview.html#)

## Module 1

  The module 2 was developed with the purpose to perform chimeric transcripts (TE-gene chimeras) detection without the need of using a reference genome and its respective annotation. This is an important approach that may be able to predict chimeric transcripts that derived from TEs that are not present in the reference genome.
  
  In the module 2, it will be performed two alignments with bowtie2 software (LANGMEAD; SALZBERG, 2012), one against transcripts and another against TE insertions, both with the parameters: -D 20 -R 3 -N 1 -L 20 -i S,1,0.50. By default, the threshold for chimeric reads to support a chimeric transcript in the module 1 is at least 2. The expression levels of genes are calculated with express (ROBERTS; PACHTER, 2013) and all low-expressed genes (FPKM < 1) for the downstream analysis. 

  In order to identify chimeric reads between TEs and gene transcripts, the alignments are converted to bed format with bedtools (QUINLAN; HALL, 2010). Then, the TE alignment output is used to create a list with all read IDs that have the mate 1 aligned against TEs, and another list with all read IDs from the mate 2. The same lists are created by using the gene alignment. All mate 1 IDs that have aligned against TEs are searched in the list of mate 2 gene IDs, which are composed by fragments that one mate has aligned against a TE and the other one against a gene transcript. In cases in which the transcript described in the reference genome does not have the TE insertion within, it will be supported only by singleton reads. Nonetheless, those cases in which the reference transcript is already described as a chimera, the chimeric reads supporting it may either come from singleton reads or concordant reads. Taking into account that the alignments are performed with end-to-end method, differentially from module 1, in the module 2 there is identification of split-reads between TEs and genes.
  
 ![Figura 1 - ChimeraTE module 1](https://i.imgur.com/YdOef5I.png)
### Inputs:

  #### 1. Stranded paired-end RNA-seq
  #### 2. Reference transcripts (.fasta)
  - In order to run ChimeraTE correctly, this fasta file **must** have a specific header pattern. All IDs must be composed firstly by the isoform ID, followed by the gene name. For instance, in _D. melanogaster_, the gene FBgn0263977 has two sequences:
  Tim17b-RA_FBgn0263977
  Tim17b-RB_FBgn0263977
  
  Note that headers "Tim17b-RA" and "Tim17b-RB" have isoform ID separated from gene name by "_". 
  This is not a usual ID format and thefore we have the script _isoform_IDs.sh_ to transform the IDs from native NCBI format to the ChimeraTE format (see details in Manual). This script may be used if you are using a genome annotation from NCBI.

  In addition, we provide here the corrected IDs for _D. melanogaster_, human (hg38), mouse (mmX) and _A. thaliana_. 

  #### 3. Reference TE insertions (.fasta)

  - This .fasta file must have only TE insertions. Be sure that they do not contains any Satellites or Low complexity repeats.
  - The .fasta file with the reference TE insertions **must** have only the TE family in the headers. For instance, if _D. melanogaster_ genome has ~4.000 DNAREP-1 TE insertions, all of them must have the header as ">DNAREP-1".

## Module 2
