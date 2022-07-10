#/bin/bash

set -e

gene_regions () {
printf "\nMeasuring TEs and genes expression\n"

awk '$3 == "gene"' "$GENE_ANNOTATION" | awk '{print $1,$4,$5,$10,$6,$7}' | sed 's/[";]//g' | sed 's/ /\t/g' > "$GENE_info"/gene_coord.bed

awk '$3 == "ncRNA"' "$GENE_ANNOTATION" | awk '{print $10}' | sed 's/[";]//g' | sort | uniq > "$GENE_info"/ncrna_ids.lst
	if [ -s "$GENE_info"/ncrna_ids.lst ]; then
		awk '$3 == "CDS"' "$GENE_ANNOTATION" | awk '{print $1,$4,$5,$10,$6,$7}' | sed 's/[";]//g' | sed 's/ /\t/g' | sort | uniq > "$GENE_info"/exon_coord1.bed
		grep -w -f "$GENE_info"/ncrna_ids.lst "$GENE_ANNOTATION" | awk '$3 == "exon"' | awk '{print $1,$4,$5,$10,$6,$7}' | sed 's/[";]//g' | sed 's/ /\t/g' | sort | uniq > "$GENE_info"/exon_coord2.bed
		cat "$GENE_info"/exon_coord1.bed "$GENE_info"/exon_coord2.bed > "$GENE_info"/exon_coord.bed
	else
		awk '$3 == "exon"' "$GENE_ANNOTATION" | awk '{print $1,$4,$5,$10,$6,$7}' | sed 's/[";]//g' | sed 's/ /\t/g' | sort | uniq > "$GENE_info"/exon_coord.bed
	fi

if [ ! -z $UTR ]; then
	echo -e "########################\n# UTR regions selected #\n########################\n"
	awk '$3 == "5UTR"' "$GENE_ANNOTATION" | awk '{print $1,$4,$5,$10,$6,$7}' | sed 's/[";]//g' | sed 's/ /\t/g' | sort | uniq > "$GENE_info"/5utr.bed
	awk '$3 == "3UTR"' "$GENE_ANNOTATION" | awk '{print $1,$4,$5,$10,$6,$7}' | sed 's/[";]//g' | sed 's/ /\t/g' | sort | uniq > "$GENE_info"/3utr.bed

	if [ -s "$GENE_info"/5utr.bed ] || [ -s "$GENE_info"/3utr.bed ]; then
		bedtools intersect -a "$ALN"/accepted_hits.bed -b "$GENE_info"/5utr.bed -wa -wb | cut -f7-12 | sort | uniq > "$GENE_info"/expressed_5utr.bed
		bedtools intersect -a "$ALN"/accepted_hits.bed -b "$GENE_info"/3utr.bed -wa -wb | cut -f7-12 | sort | uniq > "$GENE_info"/expressed_3utr.bed
	else
		echo -e "##WARNING:\nUTR option was selected, but there is no 5UTR and 3UTR annotation in the provided GTF/GFF files. Continuing without UTR analysis..."
		sleep 2
		rm "$GENE_info"/5utr.bed "$GENE_info"/3utr.bed
	fi
else
	echo "UTR regions weren't selected"
fi

if egrep -q '.gff|.gtf' <<< "$TE_ANNOTATION"; then
	awk '{print $1,$4,$5,$9,$6,$7}' "$TE_ANNOTATION" | sed 's/[";]//g' | sed 's/ /\t/g' > "$TE_info"/TE_coord.bed; else
	if grep -q '.bed'; then
	cat <<< "$TE_ANNOTATION" > "$TE_info"/TE_coord.bed; else
	echo -e "\nTE annotation formart is not valid! Exiting..." && exit 1
	fi
fi
}

expression () {
	echo "Identifying expressed genes and TEs"
	cufflinks "$PROJECT"/"$sample"/alignment/accepted_hits.bam -p "$THREADS" -G "$GENE_ANNOTATION" -o "$ALN"/fpkm/counts --quiet
	bedtools intersect -a "$ALN"/accepted_hits.bed -b "$TE_info"/TE_coord.bed -wa -wb | cut -f7-12 | sort | uniq > "$TE_info"/expressed_TEs.bed

	awk '$12 >= 1' "$ALN"/fpkm/counts/genes.fpkm_tracking | cut -f1 | sort | uniq > "$ALN"/genes_expressed_IDs.lst
	LC_ALL=C fgrep -w -f "$ALN"/genes_expressed_IDs.lst "$GENE_info"/gene_coord.bed > "$ALN"/genes_total_expressed.bed

#######	#gene_id=$(cut -f1 "$ALN"/fpkm/counts/genes.fpkm_tracking | sed '1d' | sort | uniq)
#######	#min_fpkm="1"
#######	#while -r line; do
####### #       if (( $(echo "$avg_fpkm_cov >= $min_fpkm" |bc -l) )); then
####### #       grep -w "$line" "$ALN"/fpkm_counts/results_tidy.xprs >> "$ALN"/fpkm_counts/expressed_transcripts.xprs; else
####### #      grep -w "$line" "$ALN"/fpkm_counts/results_tidy.xprs >> "$ALN"/fpkm_counts/non_exp_transcripts.xprs
#######       fi
#######	#done <<< "$gene_id"
#######	#awk '$12 >= 1' "$ALN"/fpkm/counts/genes.fpkm_tracking > "$ALN"/genes_total_expressed.bed
	bedtools intersect -a "$ALN"/accepted_hits.bed -b "$GENE_info"/gene_coord.bed -wa -wb | cut -f7-12 | sort | uniq > "$GENE_info"/expressed_genes.bed

	echo "Strand-specific expression analysis"
	bedtools intersect -a "$ALN"/fwd.bed -b "$TE_info"/TE_coord.bed -wa > "$TE_info"/TE_reads_fwd.bed
	bedtools intersect -a "$ALN"/rev.bed -b "$TE_info"/TE_coord.bed -wa > "$TE_info"/TE_reads_rev.bed
	cat "$TE_info"/TE_reads_fwd.bed "$TE_info"/TE_reads_rev.bed | cut -f4 | sed 's/[/]/*/g' | sed 's/*1//g' | sed 's/*2//g' > "$TE_info"/TE_reads.lst

	bedtools intersect -a "$ALN"/fwd.bed -b "$ALN"/genes_total_expressed.bed -wa > "$GENE_info"/gene_reads_fwd.bed
	bedtools intersect -a "$ALN"/rev.bed -b "$ALN"/genes_total_expressed.bed -wa > "$GENE_info"/gene_reads_rev.bed
	cat "$GENE_info"/gene_reads_fwd.bed "$GENE_info"/gene_reads_rev.bed | cut -f4 | sed 's/[/]/*/g' | sed 's/*1//g' | sed 's/*2//g' > "$GENE_info"/gene_reads.lst

	echo "Chimeric pairs identification"
	LC_ALL=C fgrep -F -w -f "$TE_info"/TE_reads.lst "$GENE_info"/gene_reads.lst > "$ALN"/chim_reads.lst

	LC_ALL=C fgrep -F -w -f "$ALN"/chim_reads.lst "$TE_info"/TE_reads_fwd.bed > "$TE_info"/TE_chim_readsFWD.bed
	LC_ALL=C fgrep -F -w -f "$ALN"/chim_reads.lst "$TE_info"/TE_reads_rev.bed > "$TE_info"/TE_chim_readsREV.bed

	LC_ALL=C fgrep -F -w -f "$ALN"/chim_reads.lst "$GENE_info"/gene_reads_fwd.bed > "$GENE_info"/gene_chim_readsFWD.bed
	LC_ALL=C fgrep -F -w -f "$ALN"/chim_reads.lst "$GENE_info"/gene_reads_rev.bed > "$GENE_info"/gene_chim_readsREV.bed
}

gene_regions
expression
