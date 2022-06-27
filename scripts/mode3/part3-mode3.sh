#!/bin/bash

mkdir  "$OUTPUT"/"$sample"/expression 2>/dev/null
align_and_estimate_abundance.pl --transcripts "$OUTPUT"/"$sample"/Trinity.fasta --seqType fq --left "$mate1_rep" --right "$mate2_rep" --est_method RSEM --aln_method bowtie2 --SS_lib_type "$STRANDED" --thread_count "$THREADS" --gene_trans_map  "$OUTPUT"/"$sample"/Trinity.fasta.gene_trans_map --output_dir "$OUTPUT"/"$sample"/expression --prep_reference
