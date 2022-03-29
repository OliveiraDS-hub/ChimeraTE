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
  
### Input data:

  1. Stranded paired-end RNA-seq
  2. Reference transcripts (.fasta)
  3. Reference TE insertions (.fasta)



## Module 2
