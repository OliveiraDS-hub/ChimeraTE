#!/bin/bash
set -e

term_downstream () {
while IFS= read -r line
do
	minimum=1
	scaf=$(cut -f1 <<< "$line")
	start=$(cut -f2 <<< "$line")
	end=$(cut -f3 <<< "$line")
	part=$(cut -f4,5 <<< "$line")
	fita=$(cut -f6 <<< "$line")
	if [ "$fita" == + ]; then
		num=$(($end + $WINDOW))
		echo -e "$scaf""\t""$end""\t""$num""\t""$part""\t""$fita" >> "$TMP"/genes_DOWN_window.bed
		else
		num=$(($start - $WINDOW))
			if [[ "$num" -lt "$minimum" ]]; then
                        num=1
                        fi
		echo -e "$scaf""\t""$num""\t""$start""\t""$part""\t""$fita" >> "$TMP"/genes_DOWN_window.bed
	fi
done < "$GENE_info"/expressed_genes.bed

bedtools intersect -a "$TE_info"/expressed_TEs.bed -b "$TMP"/genes_DOWN_window.bed -wb | awk '{print $10}' | sort | uniq > "$TMP"/genes_TE_DOWN.txt
}

te_terminated () {
minimum=1
grep -w -f "$TMP"/genes_TE_DOWN.txt "$GENE_info"/expressed_genes.bed > "$TMP"/gene_TE_DOWN.bed

while IFS= read -r line
do
	chr=$(cut -f1 <<< "$line")
	start=$(cut -f2 <<< "$line")
	start_reg=$(($start - "$WINDOW"))
		if [[ "$start_reg" -lt "$minimum" ]]; then
  		start_reg=1
		fi
	end=$(cut -f3 <<< "$line")
	end_reg=$(($end + "$WINDOW"))
	gene_ID=$(cut -f4 <<< "$line" | sed 's/[";]//g')
	dot=$(cut -f5 <<< "$line")
	strand=$(cut -f6 <<< "$line")
	if [ "$strand" == + ]; then
		num=$(($end + $WINDOW))
		echo -e "$chr""\t""$end""\t""$num""\t""$gene_ID""\t""$dot""\t""$strand" > "$TMP"/tmp_info.bed
		else
		num_minus=$(($start - $WINDOW))
			if [ "$num_minus" -lt "$minimum" ]; then
                        num_minus=1
                        fi
		echo -e "$chr""\t""$num_minus""\t""$start""\t""$gene_ID""\t""$dot""\t""$strand" > "$TMP"/tmp_info.bed
	fi

	bedtools intersect -a "$TE_info"/expressed_TEs.bed -b "$TMP"/tmp_info.bed -wa -wb > "$TMP"/crossing.txt
	if [ ! -s "$TMP"/crossing.txt ]; then
        continue; else
	coords_TEs=$(cut -f1-6,10 "$TMP"/crossing.txt)
	echo "$coords_TEs" > "$TMP"/TEs_at_DOWN-window.txt
		while IFS= read -r line2
		do
			chr2=$(cut -f1 <<< "$line2")
			start2=$(cut -f2 <<< "$line2")
			end_TE=$(cut -f3 <<< "$line2")
			DANIEL=$(cut -f4 <<< "$line2")
			dot2=$(cut -f5 <<< "$line2")
			strd=$(cut -f6 <<< "$line2")
			gn_id=$(cut -f7 <<< "$line2")
			if [ "$strand" == + ]; then
				if [ "$start2" -lt "$end" ]; then
				start2="$end"; fi
			else
				if [[ "$end2" -gt "$start" ]]; then
				start2="$end"; fi
			fi

  		echo -e "$chr""\t""$start""\t""$end""\t""$gene_ID""\t""$dot""\t""$strand" > "$TMP"/tmp_gene_info.bed
			echo -e "$chr2""\t""$start2""\t""$end_TE""\t""$DANIEL""\t""$dot2""\t""$strd" > "$TMP"/tmp_TE_info.bed

			if [ "$strand" == + ]; then
				reads_gene=$(bedtools intersect -a "$GENE_info"/gene_chim_readsFWD.bed -b "$TMP"/tmp_gene_info.bed)
				echo "$reads_gene" > "$TMP"/tmp_fwd_region.bed
				samt_TE=$(bedtools intersect -a "$TE_info"/TE_chim_readsFWD.bed -b "$TMP"/tmp_TE_info.bed -f "$OVERLAP" | awk '{print $4}' | sed 's/[/]/*/g' | sed 's/*1//g' | sed 's/*2//g')
				else
				reads_gene=$(bedtools intersect -a "$GENE_info"/gene_chim_readsREV.bed -b "$TMP"/tmp_gene_info.bed)
				echo "$reads_gene" > "$TMP"/tmp_rev_region.bed
				samt_TE=$(bedtools intersect -a "$TE_info"/TE_chim_readsREV.bed -b "$TMP"/tmp_TE_info.bed -f "$OVERLAP" | awk '{print $4}' | sed 's/[/]/*/g' | sed 's/*1//g' | sed 's/*2//g')
			fi

			if [ -z "$samt_TE" ]; then
  				continue; else
 					echo -e "$samt_TE" > "$TMP"/samt_TE.lst
                                        echo -e "$reads_gene" > "$TMP"/read_ids.lst
                                        counting=$(LC_ALL=C fgrep -F -w -f "$TMP"/samt_TE.lst "$TMP"/read_ids.lst | cut -f4 | sort | uniq | wc -l)
  					if [ $counting -gt 0 ]; then
  					printf "$gn_id""\t""$strand""\t""$DANIEL""\t""$strd""\t""$counting""\t""$chr":"$start"-"$end""\t""$chr2":"$start2"-"$end_TE""\n" >> "$PROJECT"/"$sample"/TE-terminated-DOWN_raw_"$sample".ct
					fi
			fi
		done < "$TMP"/TEs_at_DOWN-window.txt
	fi
done < "$TMP"/gene_TE_DOWN.bed

if [ ! -z $UTR ] ; then
	bedtools intersect -a "$TE_info"/expressed_TEs.bed -b "$GENE_info"/expressed_3utr.bed -wa -wb | cut -f7-12 | sort | uniq > "$TMP"/genes_TE_3UTR.bed
	while IFS= read -r line
	do
	echo -e "$line" > "$TMP"/tmp_info.bed
	bedtools intersect -a "$TE_info"/expressed_TEs.bed -b "$TMP"/tmp_info.bed -wa -wb > "$TMP"/crossing.txt
	if [ ! -s "$TMP"/crossing.txt ]; then
        	continue; else
       		coords_TEs=$(cut -f1,2,3,4,5,6,10 "$TMP"/crossing.txt)
      		echo "$coords_TEs" > "$TMP"/TEs_INSIDE.txt
	        chr_gene=$(cut -f1 "$TMP"/tmp_info.bed)
        	start_gene=$(cut -f2 "$TMP"/tmp_info.bed)
        	end_gene=$(cut -f3 "$TMP"/tmp_info.bed)
		strand_gene=$(cut -f6 "$TMP"/tmp_info.bed)
                	while IFS= read -r line2
                	do
                        	chr2=$(cut -f1 <<< "$line2")
                        	start2=$(cut -f2 <<< "$line2")
                        	end_TE=$(cut -f3 <<< "$line2")
                        	DANIEL=$(cut -f4 <<< "$line2")
                        	dot2=$(cut -f5 <<< "$line2")
                        	strd=$(cut -f6 <<< "$line2")
                        	gn_id=$(cut -f7 <<< "$line2")

                        echo -e "$chr2""\t""$start2""\t""$end_TE""\t""$DANIEL""\t""$dot2""\t""$strd" > "$TMP"/tmp_TE_info.bed

                        	if [ "$strd" == + ]; then
                                	reads_gene=$(bedtools intersect -a "$GENE_info"/gene_chim_readsFWD.bed -b "$TMP"/tmp_info.bed -f 0.2)
                                	samt_TE=$(bedtools intersect -a "$TE_info"/TE_chim_readsFWD.bed -b "$TMP"/tmp_TE_info.bed -f 0.2 | awk '{print $4}' | sed 's/[/]/*/g' | sed 's/*1//g' | sed 's/*2//g')
                                	else
	                                reads_gene=$(bedtools intersect -a "$GENE_info"/gene_chim_readsREV.bed -b "$TMP"/tmp_info.bed -f 0.2)
	                                samt_TE=$(bedtools intersect -a "$TE_info"/TE_chim_readsREV.bed -b "$TMP"/tmp_TE_info.bed -f 0.2 | awk '{print $4}' | sed 's/[/]/*/g' | sed 's/*1//g' | sed 's/*2//g')
        	                fi
                	        if [ -z "$samt_TE" ]; then
                        	        continue; else
                                		echo -e "$samt_TE" > "$TMP"/samt_TE.lst
                                       		echo -e "$reads_gene" > "$TMP"/read_ids.lst
                                        	counting=$(LC_ALL=C fgrep -F -w -f "$TMP"/samt_TE.lst "$TMP"/read_ids.lst | cut -f4 | sort | uniq | wc -l)
                                        	if [ $counting -gt 0 ]; then
                                        	printf "$gn_id""\t""$strand_gene""\t""$DANIEL""\t""$strd""\t""$counting""\t""$chr_gene"":""$start_gene"-"$end_gene""\t""$chr2"":""$start2"-"$end_TE""\n" >> "$PROJECT"/"$sample"/TE-terminated-3UTR_raw_"$sample".ct
                                        	fi
                        	fi
                	done < "$TMP"/TEs_INSIDE.txt
        	fi
	done < "$TMP"/genes_TE_3UTR.bed
fi
}

amb_remove-down () {
gene=$(cut -f1 "$PROJECT"/"$sample"/TE-terminated-DOWN_raw_"$sample".ct | sort | uniq)

while IFS= read -r gene_id; do
        TEs=$(grep "$gene_id" "$PROJECT"/"$sample"/TE-terminated-DOWN_raw_"$sample".ct | cut -f3 | sort | uniq)
        
        while IFS= read -r TE_id; do
                grep "$gene_id" "$PROJECT"/"$sample"/TE-terminated-DOWN_raw_"$sample".ct | grep "$TE_id" | sort -Vr -k5,5 | head -1 >> "$PROJECT"/"$sample"/TE-terminated-DOWN_"$sample".ct
        done <<< "$TEs"

done <<< "$gene"
}

amb_remove-3utr () {
gene=$(cut -f1 "$PROJECT"/"$sample"/TE-terminated-3UTR_raw_"$sample".ct | sort | uniq)

while IFS= read -r gene_id; do
        TEs=$(grep "$gene_id" "$PROJECT"/"$sample"/TE-terminated-3UTR_raw_"$sample".ct | cut -f3 | sort | uniq)

        while IFS= read -r TE_id; do
                grep "$gene_id" "$PROJECT"/"$sample"/TE-terminated-3UTR_raw_"$sample".ct | grep "$TE_id" | sort -Vr -k5,5 | head -1 >> "$PROJECT"/"$sample"/TE-terminated-3UTR_"$sample".ct
        done <<< "$TEs"

done <<< "$gene"
}

term_downstream
if [[ -s "$TMP"/genes_TE_DOWN.txt ]]; then
        genes_TE_DOWN=$(wc -l "$TMP"/genes_TE_DOWN.txt | awk '{print $1}')
        echo -e "There are $genes_TE_DOWN genes with FPKM > $FPKM with TEs at downstream region\nAnalyzing chimeric pairs for TE-terminated transcripts...\n"
        te_terminated
        if [[ -s "$PROJECT"/"$sample"/TE-terminated-DOWN_raw_"$sample".ct ]]; then
        amb_remove-down
	cp "$PROJECT"/"$sample"/TE-terminated-DOWN_"$sample".ct "$PROJECT"/tmp/; else
        echo -e "WARNING: There are no TE-terminated transcripts from TEs at downstream region in $sample"
        fi
        if [[ -s "$PROJECT"/"$sample"/TE-terminated-3UTR_raw_"$sample".ct ]]; then
        amb_remove-3utr
	cp "$PROJECT"/"$sample"/TE-terminated-3UTR_"$sample".ct "$PROJECT"/tmp/; else
        echo -e "WARNING: There are no TE-terminated transcripts from TEs at 3' UTR region in $sample"
        fi
else
        echo -e "WARNING: There are no genes with FPKM >= $FPKM and with TEs at downstream region in $sample!\n"

fi


#if [ -f "$PROJECT"/"$sample"/TE-terminated-DOWN_"$sample".ct ]; then
#        cp "$PROJECT"/"$sample"/TE-terminated-DOWN_"$sample".ct "$PROJECT"/tmp/; else
#        echo -e "There is no TE-terminated transcripts with TEs at downstream position in the "$sample" sample!"
#fi
#if [ -f "$PROJECT"/"$sample"/TE-terminated-3UTR_"$sample".ct ]; then
#        cp "$PROJECT"/"$sample"/TE-terminated-3UTR_"$sample".ct "$PROJECT"/tmp/; else
#        echo -e "There is no TE-terminated transcripts with TEs at 3' UTR in the "$sample" sample!"
#fi

rm "$TMP"/genes_DOWN_window.bed; rm "$PROJECT"/"$sample"/*_raw_*
