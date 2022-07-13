#!/bin/bash

all_chim_genes=$(cut -f1 "$READS"/genes_chim.bed | sort | uniq)

echo "Capturing sequence from chimeric read pairs"

seqtk subseq "$mate1" "$READS"/all_chim_IDs.lst > "$READS"/reads_R1.fq
seqtk subseq "$mate2" "$READS"/all_chim_IDs.lst > "$READS"/reads_R2.fq

seqtk subseq "$mate1" "$READS"/all_chim_IDs_1.lst > "$READS"/reads_TE_R1.fq
seqtk subseq "$mate2" "$READS"/all_chim_IDs_2.lst > "$READS"/reads_TE_R2.fq

echo "Recovering chimeric transcripts"

while IFS= read -r gene_id
do
  gene_info=$(awk -v id="$gene_id" '$1 == id' "$READS"/genes_chim.bed)
	read_gene_ids=$(cut -f4 <<< "$gene_info" | sort | uniq)
#	#read_gene=$(cut -f4 <<< "$gene_info" | head -1)
#	#read_TE=$(LC_ALL=C fgrep -w "$read_gene_ids" "$READS"/TEs_chim.bed | cut -f4 | head -1)
	check=$(LC_ALL=C fgrep -w "$read_gene_ids" "$READS"/TEs_chim.bed | cut -f1 | sort | uniq -c | sort -V | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | head -1 | awk '{print $2,$1}'| sed 's/ /\t/g')
    if [ ! -z "$check" ]; then
      TE_ID=$(LC_ALL=C fgrep -w "$read_gene_ids" "$READS"/TEs_chim.bed | cut -f1 | sort | uniq -c | sort -V | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | sort -r -k1,1 | head -1 | cut -f2)
			chim_cov=$(LC_ALL=C fgrep -w "$read_gene_ids" "$READS"/TEs_chim.bed | cut -f1 | sort | uniq -c | sort -V | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | sort -r -k1,1 | head -1 | cut -f1)
      if [ ! -z "$chim_cov" ]; then
        printf "$gene_id""\t""$TE_ID""\t""$chim_cov""\n" >> "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct
      fi
    fi
done <<< "$all_chim_genes"

					#TE_ID_seqtk=$(grep -w "$read_TE" "$READS"/all_chim_IDs_1.lst)
						#if [ ! -z "$TE_ID_seqtk" ]; then
						#	echo -e "$TE_ID_seqtk" > "$READS"/read_TE.lst
						#	TE_seq=$(seqtk subseq "$READS"/reads_TE_R1.fq "$READS"/read_TE.lst | head -2 | tail -1)
						#	TE_read_strd="-"
						#	gene_seq=$(seqtk subseq "$READS"/reads_R2.fq "$READS"/read_TE.lst | head -2 | tail -1)
						#	gene_read_strd="+"
						#else
						#	TE_ID_seqtk=$(grep -w "$read_TE" "$READS"/all_chim_IDs_2.lst)
						#	echo -e "$TE_ID_seqtk" > "$READS"/read_TE.lst
						#	TE_seq=$(seqtk subseq "$READS"/reads_TE_R2.fq "$READS"/read_TE.lst | head -2 | tail -1)
						#	TE_read_strd="+"
                                                 #       gene_seq=$(seqtk subseq "$READS"/reads_R1.fq "$READS"/read_TE.lst | head -2 | tail -1)
						#	gene_read_strd="-"
						#fi

genes=$(cut -f1 "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct | sed 's/_/\t/g' | cut -f2 | sort | uniq)
while IFS= read -r line; do
	TE_ID=$(grep "$line" "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct | cut -f2 | head -1)
	total_cov=$(grep "$line" "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct | awk '{sum+=$3;}END{print sum;}')

#	if [ "$total_cov" -gt "$CUTOFF" ]; then
		isoforms=$(grep "$line" "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct | cut -f1 | sed -z 's/\n/,/g;s/,$/\n/')
		printf "$line""\t""$TE_ID""\t""$total_cov""\t""$isoforms""\n" >> "$PROJECT"/"$sample"/tmp_dir/"$sample"_chimTEs_final.ct
#	fi
done <<< "$genes"


#min_co_exp="10"
#min_fpkm="1"

#printf "\nCalculating chimeric reads coverage... \n"

#sort -k1,1 "$sample"/chimTEs_final.ct | cut -f1-3 > "$ALN"/chimTEs_sorted.ct
#IDs=$(cut -f1 "$ALN"/chimTEs_sorted.ct)
#grep -w "$IDs" "$ALN"/fpkm_counts/expressed_transcripts.xprs | awk '{print $2,$5,$13}' | awk '{$3 = sprintf("%.3f", $3)} 1' | sed 's/ /\t/g' | sort -k1,1 | cut -f2,3 > "$ALN"/fpkm_counts/fpkm_chimeras.xprs
#chim_expression=$(paste "$ALN"/chimTEs_sorted.ct "$ALN"/fpkm_counts/fpkm_chimeras.xprs)

#while IFS= read -r line; do
#	fpkm=$(cut -f5 <<< "$line")
#	DONT USE IT		#if (( $(echo "$fpkm >= $min_fpkm" |bc -l) )); then
#			cov_chim=$(cut -f3 <<< "$line")
#			cov_genes=$(cut -f4 <<< "$line")

#			co_exp=$(expr "$fpkm"*"$cov_chim"/"$cov_genes" |bc -l )
#			co_exp=$(awk -v fpkm=$fpkm -v cov_chim=$cov_chim -v cov_genes=$cov_genes 'BEGIN { print  ( fpkm * cov_chim / cov_genes ) }')
#			printf "$line""\t""$co_exp""\n" >> "$sample"/chimTE_fpkm_results.ct

#done <<< "$chim_expression"

#printf "\nChimeric reads coverage done! \n"

#cut -f2 "$proj"/chimTEs_final.ct | sort | uniq -c | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | awk '{print $2,$1}' | sed 's/ /\t/g' > "$proj"/plot/TE-freq.tsv

#cp scripts/plotting/plot.R "$proj"/plot/

#cd "$proj"/plot/
#Rscript plot.R





#echo "REPLICABILITY"






#echo "Process successfully finished!"
