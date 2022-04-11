#!/bin/bash

set -e

alignments () {
	printf "\nRunning TEs alignment...\n"

	if [ -f indexes/tes_index_"$proj".1.bt2 ]; then
		printf "\nBowtie2 index for $proj has been found\n"; else
		bowtie2-build "$input_TEs" indexes/tes_index_"$proj" --quiet && printf "\n$clock\n" && printf "\nBowtie2 index created with $input_TEs\n"
	fi

	bowtie2 -x indexes/tes_index_"$proj" -1 "$mate1" -2 "$mate2" -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p 24 -S "$ALN"/tes.sam

	printf "\nTEs alignment is done! \n"

	printf "\nRunning genes alignment...\n"

	if [ -f indexes/genes_index_"$proj".1.bt2 ]; then
		printf "\nBowtie2 index for $proj has been found\n"; else
		bowtie2-build "$input_genes" indexes/genes_index_"$proj" --quiet && printf "\n$clock\n" && printf "\nBowtie2 index created with $input_genes\n"
	fi
        
	bowtie2 -x indexes/genes_index_"$proj" -1 "$mate1" -2 "$mate2" -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p 24 -S "$ALN"/genes.sam

	printf "\nGenes alignment is done! \n"
}

conversion () {
	for SAMF in "$ALN"/*.sam; do
	name="${SAMF%.*}"
	output="${name}"
	samtools view -bS "$SAMF" > "$output".bam
	bedtools bamtobed -i "$output".bam > "$output"_total.bed
	done
	printf "\n******Samtools has finished all conversions!******"
}

expression () {
	printf "\nCalculating FPKM to the transcripts..."
	express -o "$ALN"/fpkm_counts -O 1 --output-align-prob --no-bias-correct --"$strandness" "$input_genes" "$ALN"/genes.bam
	printf "\nFPKM done!"

	printf "\nRemoving transcripts with FPKM lower than 1..."

	min_fpkm="1"
	awk '{$13 = sprintf("%.3f", $13)} 1' "$ALN"/fpkm_counts/results.xprs | sed 's/[ _]/\t/g' > "$ALN"/fpkm_counts/results_tidy.xprs
	gene_id=$(cut -f2 "$ALN"/fpkm_counts/results.xprs | sed '1d' | sed 's/_/\t/g' | cut -f2 | sort | uniq)
	
	while read -r line;do
		avg_fpkm_cov=$(grep -w "$line" "$ALN"/fpkm_counts/results_tidy.xprs | awk '{print $14}' | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }')

		if (( $(echo "$avg_fpkm_cov >= $min_fpkm" |bc -l) )); then
		grep -w "$line" "$ALN"/fpkm_counts/results_tidy.xprs >> "$ALN"/fpkm_counts/expressed_transcripts.xprs; else
		grep -w "$line" "$ALN"/fpkm_counts/results_tidy.xprs >> "$ALN"/fpkm_counts/non_exp_transcripts.xprs
		fi
	done <<< "$gene_id"
	
	cut -f3 "$ALN"/fpkm_counts/non_exp_transcripts.xprs | sort | uniq > "$ALN"/fpkm_counts/non_exp_transcripts.lst
	sed 's/_/ /g' "$ALN"/genes_total.bed | LC_ALL=C fgrep -v -w -f "$ALN"/fpkm_counts/non_exp_transcripts.lst | sed 's/ /_/g' > "$ALN"/genes_total_expressed.bed

}

chim_IDs () {
	printf "\nIdentifying chimeric read IDs... \n"
	cut -f4 "$ALN"/tes_total.bed | grep '/1' | sed 's|/1||' > "$ALN"/TEs_R1_IDs.lst
	cut -f4 "$ALN"/tes_total.bed | grep '/2' | sed 's|/2||' > "$ALN"/TEs_R2_IDs.lst
	cut -f4 "$ALN"/genes_total_expressed.bed | grep '/1' | sed 's|/1||' > "$ALN"/genes_R1_IDs.lst
	cut -f4 "$ALN"/genes_total_expressed.bed | grep '/2' | sed 's|/2||' > "$ALN"/genes_R2_IDs.lst

	LC_ALL=C fgrep -w -f "$ALN"/TEs_R1_IDs.lst "$ALN"/genes_R2_IDs.lst | sort | uniq > "$READS"/all_chim_IDs_1.lst
	LC_ALL=C fgrep -w -f "$ALN"/TEs_R2_IDs.lst "$ALN"/genes_R1_IDs.lst | sort | uniq > "$READS"/all_chim_IDs_2.lst
	cat "$READS"/all_chim_IDs_1.lst "$READS"/all_chim_IDs_2.lst | sort | uniq > "$READS"/all_chim_IDs.lst

	LC_ALL=C fgrep -w -f "$READS"/all_chim_IDs.lst "$ALN"/tes_total.bed | sed 's|/2||; s|/1||' > "$READS"/TEs_chim.bed
	LC_ALL=C fgrep -w -f "$READS"/all_chim_IDs.lst "$ALN"/genes_total_expressed.bed | sed 's|/2||; s|/1||' > "$READS"/genes_chim.bed

	printf "\ndone! \n"
}

alignments
conversion
expression
chim_IDs
