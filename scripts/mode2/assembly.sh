#!/bin/bash

#set -e

assembly_trin () {
#assembly de novo transcripts
if [[ $STRANDED == "rf-stranded" ]]; then
  STRANDED="RF"; else
  if [[ $STRANDED == "fr-stranded" ]]; then
    STRANDED="FR"
  fi
fi
Trinity --seqType fq --SS_lib_type "$STRANDED" --max_memory "$RAM"G --left "$mate1" --right "$mate2" --CPU "$THREADS" --output "$TRINITY_OUT" 1> /dev/null
echo -e "${GREEN}DONE!!${NC}\n"; sleep 2
rm -R "$TRINITY_OUT"/chrysalis "$TRINITY_OUT"/insilico_read_normalization "$TRINITY_OUT"/read_partitions/ "$TRINITY_OUT"/scaffolding_entries.sam "$TRINITY_OUT"/both* #"$TRINITY_OUT"/jellyfish.kmers.fa
}

massive_data () {
#alignment
echo -ne "Creating bowtie2-index for assembled transcripts..."
bowtie2-build "$TRINITY_OUT"/Trinity.fasta "$TRINITY_OUT"/trinity-assembly --threads "$THREADS" >> /dev/null 2>&1; echo -e "${GREEN}DONE!${NC}\n"

echo -e "Performing bowtie2 alignment..."
bowtie2 -x "$TRINITY_OUT"/trinity-assembly -1 "$MATE1" -2 "$MATE2" -S "$TRINITY_OUT"/accepted_hits.sam -p "$THREADS"; echo -e "${GREEN}DONE!${NC}\n"

#conversion
echo -e "Converting sam to bed..."
samtools view -@ "$THREADS" -bS "$TRINITY_OUT"/accepted_hits.sam > "$TRINITY_OUT"/accepted_hits.bam; rm "$TRINITY_OUT"/accepted_hits.sam
bedtools bamtobed -i "$TRINITY_OUT"/accepted_hits.bam > "$TRINITY_OUT"/accepted_hits.bed; rm "$TRINITY_OUT"/accepted_hits.bam; echo -e "${GREEN}DONE!${NC}\n"

#masking transcripts
if egrep -q '.fasta|.fa' <<< "$REF_TEs"; then
	echo -ne "Running RepeatMasker with built-in library $REF_TEs..."
	RepeatMasker "$TRINITY_OUT"/Trinity.fasta -lib "$REF_TEs" -cutoff 225 -nolow -norna -s -par "$THREADS" >> /dev/null 2>&1; echo -e "${GREEN}DONE!${NC}\n"; else
	echo -ne "Running RepeatMasker with species = $REF_TEs..."
	RepeatMasker "$TRINITY_OUT"/Trinity.fasta -species "$REF_TEs" -cutoff 225 -nolow -norna -s -par "$THREADS" >> /dev/null 2>&1; echo -e "${GREEN}DONE!${NC}\n"
fi
}

############################################

managing_files () {

echo -ne "Identifying chimeric reads..."
cat "$TRINITY_OUT"/Trinity.fasta.out | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | tail -n +4 | awk 'BEGIN {FS=OFS="\t"} NR >=1 {print $0, $7 - $6}' | awk -v TE_length=$TE_length '$16 >= TE_length' > "$TRINITY_OUT"/trinity.out.tsv
awk '{print $5,$6,$7,$10,$16,$9}' "$TRINITY_OUT"/trinity.out.tsv | sed 's/ /\t/g' > "$TRINITY_OUT"/trinity.out.bed

#chimeric reads
bedtools intersect -a "$TRINITY_OUT"/accepted_hits.bed -b "$TRINITY_OUT"/trinity.out.bed -f "$OVERLAP" -wa -wb > "$TRINITY_OUT"/overlapped.tsv

cut -f4 "$TRINITY_OUT"/overlapped.tsv | sed 's/..$//g' | sort | uniq -c | awk '$1 == "1"' | awk '{print $2}' > "$TRINITY_OUT"/singletons_TEs1.lst
cut -f4 "$TRINITY_OUT"/overlapped.tsv > "$TRINITY_OUT"/all_reads.lst

grep -v -w -f "$TRINITY_OUT"/singletons_TEs1.lst "$TRINITY_OUT"/all_reads.lst > "$TRINITY_OUT"/putative_reads.lst
sed 's/..$//g' "$TRINITY_OUT"/putative_reads.lst | sort | uniq > "$TRINITY_OUT"/pairing.lst

while IFS= read -r IDs; do
	pair=$(grep -w "$IDs" "$TRINITY_OUT"/putative_reads.lst)
	mate1=$(grep '/1' <<< "$pair")
	mate2=$(grep '/2' <<< "$pair")
	if [[ ! -z $mate1 && ! -z $mate2 ]]; then
		continue ; else
		echo "$IDs" >> "$TRINITY_OUT"/1mate_2TEs.lst
	fi
done < "$TRINITY_OUT"/pairing.lst

if [[ -f "$TRINITY_OUT"/1mate_2TEs.lst ]]; then
	cat "$TRINITY_OUT"/1mate_2TEs.lst "$TRINITY_OUT"/singletons_TEs1.lst > "$TRINITY_OUT"/chim_reads_search.lst; else
	mv "$TRINITY_OUT"/singletons_TEs1.lst "$TRINITY_OUT"/chim_reads_search.lst
fi

if [[ ! -f "$TRINITY_OUT"/chim_reads_search.lst ]]; then
	echo -e "\nNO CHIMERIC READS IDENTIFIED!\nFinishing analysis for "$sample"" && exit 1
fi

grep -w -f "$TRINITY_OUT"/chim_reads_search.lst "$TRINITY_OUT"/overlapped.tsv | cut -f1,4,7,10 > "$TRINITY_OUT"/all_candidates.tsv

grep -w -f "$TRINITY_OUT"/chim_reads_search.lst "$TRINITY_OUT"/accepted_hits.bed > "$TRINITY_OUT"/transcripts_chim.lst
echo -e "\t${GREEN}===> DONE!${NC}\n"
}

assembly_trin
massive_data
managing_files
