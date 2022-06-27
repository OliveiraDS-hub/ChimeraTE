#!/bin/bash

cut -f1 "$PROJECT"/TE-exonized.ct | sort | uniq > "$PROJECT"/TE-exonized_genes.lst

while IFS= read -r line
do
	TEs_number=$(grep -F -w "$line" "$PROJECT"/TE-exonized.ct | cut -f3 | sort | uniq | wc -l)
	TE_strd=$(grep -F -w "$line" "$PROJECT"/TE-exonized.ct | cut -f4 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
	TE_IDs=$(grep -F -w "$line" "$PROJECT"/TE-exonized.ct | cut -f3 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
	TE_pos=$(grep -F -w "$line" "$PROJECT"/TE-exonized.ct | cut -f7 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
	EXON_number=$(grep -F -w "$line" "$PROJECT"/TE-exonized.ct | wc -l)

	GENE_ID=$(grep -F -w "$line" "$PROJECT"/TE-exonized.ct | cut -f1 | head -1)
	GENE_strd=$(grep -F -w "$line" "$PROJECT"/TE-exonized.ct | cut -f2 | head -1)
	cov=$(grep -F -w "$line" "$PROJECT"/TE-exonized.ct | cut -f5)
	if [ "$EXON_number" -gt 1 ]; then
		avg=$(awk '{ total += $1 } END { print total/NR }' <<< "$cov")
		exons=$(grep -F -w "$GENE_ID" "$PROJECT"/TE-exonized.ct | cut -f6 | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
		printf "$GENE_ID""\t""$GENE_strd""\t""$TEs_number""\t""$TE_IDs""\t""$TE_strd""\t""$TE_pos""\t""$EXON_number""\t""$exons""\t""$avg""\n" >> "$PROJECT"/TE-exonized_merging.ct; else
		exons=$(grep -F -w "$GENE_ID" "$PROJECT"/TE-exonized.ct | cut -f6)
		printf "$GENE_ID""\t""$GENE_strd""\t""$TEs_number""\t""$TE_IDs""\t""$TE_strd""\t""$TE_pos""\t""$EXON_number""\t""$exons""\t""$cov""\n" >> "$PROJECT"/TE-exonized_merging.ct
	fi
done < "$PROJECT"/TE-exonized_genes.lst

cutoff="2"

while IFS= read -r line
do
	cov=$(cut -f7 <<< "$line")
	if (( $(bc <<< "$cov >= $cutoff" ) )); then
		echo "$line" >> "$PROJECT"/TE-exonized_final.ct
	fi
done < "$PROJECT"/TE-exonized_merging.ct

sed -i '1i\Gene_ID\tTE_number\tTE_IDs\tTE_position\tExon_number\tExon_position\tChimeric pairs coverage' "$PROJECT"/TE-exonized_final.ct
