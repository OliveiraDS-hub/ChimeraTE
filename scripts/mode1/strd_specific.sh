#!/bin/bash

set -e
GENOME_file=$(basename "$GENOME" | sed 's/.fasta//g' | sed 's/.fa//g')

alignment () {
if [[ ! -d "$PROJECT"/indexes/ ]]; then
	mkdir "$PROJECT"/indexes/ 2>/dev/null
fi

if [ -f "$PROJECT"/indexes/"$GENOME_file"_index.1.bt2 ]; then
	echo -e "\nBowtie2 index for $GENOME_file has been found"; else
	echo -e "Creating Bowtie2 index..."
	bowtie2-build "$GENOME" "$PROJECT"/indexes/"$GENOME_file"_index --threads "$THREADS" >> /dev/null 2>&1 && echo -e "\nBowtie2 index created with $GENOME_file"; echo -e "${GREEN}DONE!!${NC}\n"
fi

	echo -e "Bowtie2 alignment...\n"
	bowtie2 -x "$PROJECT"/indexes/"$GENOME_file"_index -1 "$mate1_rep" -2 "$mate2_rep" -S "$ALN"/accepted_hits.sam --threads "$THREADS"

	echo -ne "\nConverting sam to bed...\n"
	samtools view -@ "$THREADS" -bS "$ALN"/accepted_hits.sam > "$ALN"/accepted_hits.bam; rm "$ALN"/accepted_hits.sam
	samtools sort -@ "$THREADS" -o "$ALN"/accepted_hits_sorted.bam "$ALN"/accepted_hits.bam; rm "$ALN"/accepted_hits.bam
	mv "$ALN"/accepted_hits_sorted.bam "$ALN"/accepted_hits.bam

	samtools index "$ALN"/accepted_hits.bam
	bedtools bamtobed -i "$ALN"/accepted_hits.bam > "$ALN"/accepted_hits.bed
	BAM_FILE="$ALN"/accepted_hits.bam
	echo -e "${GREEN}DONE!!${NC}\n"
}

management () {

	echo -ne "Identifying strand-specific reads...\n"
	# Forward strand
	samtools view -@ "$THREADS" -b -f 128 -F 16 "$BAM_FILE" > "$PROJECT"/"$sample"/alignment/fwd1_f.bam
	samtools index "$PROJECT"/"$sample"/alignment/fwd1_f.bam
	samtools view -@ "$THREADS" -b -f 80 "$BAM_FILE" > "$PROJECT"/"$sample"/alignment/fwd2_f.bam
	samtools index "$PROJECT"/"$sample"/alignment/fwd2_f.bam

	# Combine alignments that originate on the forward strand.
	samtools merge -@ "$THREADS" -f "$PROJECT"/"$sample"/alignment/fwd.bam "$PROJECT"/"$sample"/alignment/fwd1_f.bam "$PROJECT"/"$sample"/alignment/fwd2_f.bam
	rm "$PROJECT"/"$sample"/alignment/*f.bam*
	samtools index "$PROJECT"/"$sample"/alignment/fwd.bam

	# Reverse strand
	samtools view -@ "$THREADS" -b -f 144 "$BAM_FILE" > "$PROJECT"/"$sample"/alignment/rev1_r.bam
	samtools index "$PROJECT"/"$sample"/alignment/rev1_r.bam
	samtools view -@ "$THREADS" -b -f 64 -F 16 "$BAM_FILE" > "$PROJECT"/"$sample"/alignment/rev2_r.bam
	samtools index "$PROJECT"/"$sample"/alignment/rev2_r.bam

	# Combine alignments that originate on the reverse strand.
	samtools merge -@ "$THREADS" -f "$PROJECT"/"$sample"/alignment/rev.bam "$PROJECT"/"$sample"/alignment/rev1_r.bam "$PROJECT"/"$sample"/alignment/rev2_r.bam
	rm "$PROJECT"/"$sample"/alignment/*r.bam*
	samtools index "$PROJECT"/"$sample"/alignment/rev.bam

	bedtools bamtobed -i "$ALN"/fwd.bam > "$ALN"/fwd.bed
	bedtools bamtobed -i "$ALN"/rev.bam > "$ALN"/rev.bed
	echo -e "${GREEN}DONE!!${NC}"
}
alignment
management
