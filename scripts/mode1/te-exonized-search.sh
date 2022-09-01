#!/bin/bash

set -e

init_exon () {
bedtools intersect -a "$TE_info"/expressed_TEs.bed -b "$GENE_info"/expressed_genes.bed -wb -nonamecheck | awk '{print $10}' | sort | uniq > "$TMP"/genes_TE_INSIDE.txt

grep -w -f "$TMP"/genes_TE_INSIDE.txt "$GENE_info"/exon_coord.bed > "$TMP"/gene_TE_INSIDE.bed
}

te_exonized () {
while IFS= read -r line
do
	echo -e "$line" > "$TMP"/tmp_info.bed
	bedtools intersect -a "$TE_info"/expressed_TEs.bed -b "$TMP"/tmp_info.bed -wa -wb -nonamecheck > "$TMP"/crossing.txt
	if [ ! -s "$TMP"/crossing.txt ]; then
        continue; else
	coords_TEs=$(cut -f1-6,10 "$TMP"/crossing.txt)
	echo "$coords_TEs" > "$TMP"/TEs_INSIDE.txt
	chr_gene=$(cut -f1 "$TMP"/tmp_info.bed)
        start_gene=$(cut -f2 "$TMP"/tmp_info.bed)
        end_gene=$(cut -f3 "$TMP"/tmp_info.bed)
	strand=$(cut -f6 <<< "$line")
		while IFS= read -r line2
		do
			chr2=$(cut -f1 <<< "$line2")
			start2=$(cut -f2 <<< "$line2")
			end_TE=$(cut -f3 <<< "$line2")
			TE_ID=$(cut -f4 <<< "$line2")
			dot2=$(cut -f5 <<< "$line2")
			strd=$(cut -f6 <<< "$line2")
			gn_id=$(cut -f7 <<< "$line2")

			echo -e "$chr2""\t""$start2""\t""$end_TE""\t""$TE_ID""\t""$dot2""\t""$strd" > "$TMP"/tmp_TE_info.bed

			if [ "$strand" == + ]; then
				reads_gene=$(bedtools intersect -a "$GENE_info"/gene_chim_readsFWD.bed -b "$TMP"/tmp_info.bed -nonamecheck)
				echo "$reads_gene" > "$TMP"/tmp_fwd_region.bed
				samt_TE=$(bedtools intersect -a "$TE_info"/TE_chim_readsFWD.bed -b "$TMP"/tmp_TE_info.bed -f "$OVERLAP" -nonamecheck | awk '{print $4}' | sed 's/[/]/*/g' | sed 's/*1//g' | sed 's/*2//g')
				else
				reads_gene=$(bedtools intersect -a "$GENE_info"/gene_chim_readsREV.bed -b "$TMP"/tmp_info.bed -nonamecheck)
				echo "$reads_gene" > "$TMP"/tmp_rev_region.bed
				samt_TE=$(bedtools intersect -a "$TE_info"/TE_chim_readsREV.bed -b "$TMP"/tmp_TE_info.bed -f "$OVERLAP" -nonamecheck | awk '{print $4}' | sed 's/[/]/*/g' | sed 's/*1//g' | sed 's/*2//g')
			fi
			if [ -z "$samt_TE" ]; then
  				continue; else
 					echo -e "$samt_TE" > "$TMP"/samt_TE.lst
                                        echo -e "$reads_gene" > "$TMP"/read_ids.lst
                                        counting=$(LC_ALL=C fgrep -F -w -f "$TMP"/samt_TE.lst "$TMP"/read_ids.lst | cut -f4 | sort | uniq | wc -l)
  					if [ $counting -gt 0 ]; then
  					printf "$gn_id""\t""$strand""\t""$TE_ID""\t""$strd""\t""$counting""\t""$chr_gene"":""$start_gene"-"$end_gene""\t""$chr2"":""$start2"-"$end_TE""\n" >> "$PROJECT"/"$sample"/TE-exonized_raw_"$sample".ct
					fi
			fi

		done < "$TMP"/TEs_INSIDE.txt
	fi
done < "$TMP"/gene_TE_INSIDE.bed
}

amb_remove () {
gene=$(cut -f1 "$PROJECT"/"$sample"/TE-exonized_raw_"$sample".ct | sort | uniq)

while IFS= read -r gene_id; do
	TEs=$(grep "$gene_id" "$PROJECT"/"$sample"/TE-exonized_raw_"$sample".ct | cut -f3 | sort | uniq)

	while IFS= read -r TE_id; do
		grep "$gene_id" "$PROJECT"/"$sample"/TE-exonized_raw_"$sample".ct | grep "$TE_id" | sort -Vr -k5,5 | head -1 >> "$PROJECT"/"$sample"/TE-exonized_"$sample".ct
	done <<< "$TEs"

done <<< "$gene"
}

init_exon
if [[ -s "$TMP"/genes_TE_INSIDE.txt ]]; then
        genes_TE_DOWN=$(wc -l "$TMP"/genes_TE_INSIDE.txt | awk '{print $1}')
        echo -ne "There are $genes_TE_DOWN genes with FPKM > $FPKM that have TEs inside genes\nAnalyzing chimeric pairs for TE-exonized transcripts...\n"
        te_exonized
				echo -e "${GREEN}DONE!!${NC}"
        if [[ -s "$PROJECT"/"$sample"/TE-exonized_raw_"$sample".ct ]]; then
        amb_remove
        cp "$PROJECT"/"$sample"/TE-exonized_"$sample".ct "$PROJECT"/tmp/; else
        echo -e "WARNING: There are no TE-exonized transcripts in $sample"
        fi
else
        echo -e "WARNING: There are no genes with FPKM >= $FPKM with TEs inside it in $sample!\n"
fi

rm "$PROJECT"/"$sample"/*_raw_*
echo -e "\n#################################################\n# ChimeraTE analysis with $sample has finished! #\n#################################################\n"
