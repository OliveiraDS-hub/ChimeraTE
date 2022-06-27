#!/bin/bash

set -e

#massive_data () {
#assembly de novo transcripts
#Trinity --seqType fq --SS_lib_type "$STRANDED" --max_memory "$RAM"G --left "$mate1_rep" --right "$mate2_rep" --CPU "$THREADS" #--output "$OUTPUT"/"$sample"
#mv trinity_out_dir "$OUTPUT"/"$sample"
#cp "$OUTPUT"/"$sample"/trinity_out_dir/Trinity.fasta "$OUTPUT"/"$sample"/trinity_out_dir/Trinity.fasta.gene_trans_map "$OUTPUT"/"$sample"
#alignment
#bowtie2-build "$OUTPUT"/"$sample"/Trinity.fasta "$OUTPUT"/"$sample"/trinity-assembly --quiet --threads "$THREADS"
#bowtie2 -x "$OUTPUT"/"$sample"/trinity-assembly -1 "$MATE1" -2 "$MATE2" -S "$OUTPUT"/"$sample"/accepted_hits.sam -p "$THREADS"

#conversion
#samtools view -@ "$THREADS" -bS "$OUTPUT"/"$sample"/accepted_hits.sam > "$OUTPUT"/"$sample"/accepted_hits.bam
#bedtools bamtobed -i "$OUTPUT"/"$sample"/accepted_hits.bam > "$OUTPUT"/"$sample"/accepted_hits.bed

#masking transcripts
if egrep -q '.fasta|.fa' <<< "$REF_TEs"; then
	RepeatMasker "$OUTPUT"/"$sample"/Trinity.fasta -library "$REF_TEs" -cutoff 225 -nolow -norna -s -par "$THREADS"; else
	RepeatMasker "$OUTPUT"/"$sample"/Trinity.fasta -species "$REF_TEs" -cutoff 225 -nolow -norna -s -par "$THREADS"
fi
#}

############################################

#managing_files () {
cat "$OUTPUT"/"$sample"/Trinity.fasta.out | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | tail -n +4 | awk 'BEGIN {FS=OFS="\t"} NR >=1 {print $0, $7 - $6}' | awk -v TE_length=$TE_length '$16 >= TE_length' > "$OUTPUT"/"$sample"/trinity.out.tsv
awk '{print $5,$6,$7,$10,$16,$9}' "$OUTPUT"/"$sample"/trinity.out.tsv | sed 's/ /\t/g' > "$OUTPUT"/"$sample"/trinity.out.bed

#chimeric reads
bedtools intersect -a "$OUTPUT"/"$sample"/accepted_hits.bed -b "$OUTPUT"/"$sample"/trinity.out.bed -f "$OVERLAP" -wa -wb > "$OUTPUT"/"$sample"/overlapped.tsv

cut -f4 "$OUTPUT"/"$sample"/overlapped.tsv | sed 's/..$//g' | sort | uniq -c | awk '$1 == "1"' | awk '{print $2}' > "$OUTPUT"/"$sample"/singletons_TEs1.lst
cut -f4 "$OUTPUT"/"$sample"/overlapped.tsv > "$OUTPUT"/"$sample"/all_reads.lst

grep -v -w -f "$OUTPUT"/"$sample"/singletons_TEs1.lst "$OUTPUT"/"$sample"/all_reads.lst > "$OUTPUT"/"$sample"/putative_reads.lst
sed 's/..$//g' "$OUTPUT"/"$sample"/putative_reads.lst | sort | uniq > "$OUTPUT"/"$sample"/pairing.lst

echo "finding all reads..."
while IFS= read -r IDs; do
	pair=$(grep -w "$IDs" "$OUTPUT"/"$sample"/putative_reads.lst)
	mate1=$(grep '/1' <<< "$pair")
	mate2=$(grep '/2' <<< "$pair")
	if [[ ! -z $mate1 && ! -z $mate2 ]]; then
		continue ; else
		echo "$IDs" >> "$OUTPUT"/"$sample"/1mate_2TEs.lst
	fi
done < "$OUTPUT"/"$sample"/pairing.lst

if [[ -f "$OUTPUT"/"$sample"/1mate_2TEs.lst ]]; then
	cat "$OUTPUT"/"$sample"/1mate_2TEs.lst "$OUTPUT"/"$sample"/singletons_TEs1.lst > "$OUTPUT"/"$sample"/chim_reads_search.lst; else
	mv "$OUTPUT"/"$sample"/singletons_TEs1.lst "$OUTPUT"/"$sample"/chim_reads_search.lst
fi

if [[ ! -f "$OUTPUT"/"$sample"/chim_reads_search.lst ]]; then
	echo -e "\nNO CHIMERIC READS IDENTIFIED!\nFinishing analysis for "$sample"" && exit 1
fi

grep -w -f "$OUTPUT"/"$sample"/chim_reads_search.lst "$OUTPUT"/"$sample"/overlapped.tsv | cut -f1,4,7,10 > "$OUTPUT"/"$sample"/all_candidates.tsv

grep -w -f "$OUTPUT"/"$sample"/chim_reads_search.lst "$OUTPUT"/"$sample"/accepted_hits.bed > "$OUTPUT"/"$sample"/transcripts_chim.lst
#}

#massive_data
#managing_files
