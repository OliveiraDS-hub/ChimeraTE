#!/bin/bash

all_chim_genes=$(cut -f1 "$READS"/genes_chim.bed | sort | uniq)

echo -e "Capturing sequence from chimeric read pairs..."
seqtk subseq "$mate1" "$READS"/all_chim_IDs.lst > "$READS"/reads_R1.fq
seqtk subseq "$mate2" "$READS"/all_chim_IDs.lst > "$READS"/reads_R2.fq
seqtk subseq "$mate1" "$READS"/all_chim_IDs_1.lst > "$READS"/reads_TE_R1.fq
seqtk subseq "$mate2" "$READS"/all_chim_IDs_2.lst > "$READS"/reads_TE_R2.fq

while IFS= read -r gene_id
do
  gene_info=$(awk -v id="$gene_id" '$1 == id' "$READS"/genes_chim.bed)
	read_gene_ids=$(cut -f4 <<< "$gene_info" | sort | uniq)
	check=$(LC_ALL=C fgrep -w "$read_gene_ids" "$READS"/TEs_chim.bed | cut -f1 | sort | uniq -c | sort -V | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | head -1 | awk '{print $2,$1}'| sed 's/ /\t/g')
    if [ ! -z "$check" ]; then
      TE_ID=$(LC_ALL=C fgrep -w "$read_gene_ids" "$READS"/TEs_chim.bed | cut -f1 | sort | uniq -c | sort -V | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | sort -r -k1,1 | head -1 | cut -f2)
			chim_cov=$(LC_ALL=C fgrep -w "$read_gene_ids" "$READS"/TEs_chim.bed | cut -f1 | sort | uniq -c | sort -V | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | sort -r -k1,1 | head -1 | cut -f1)
      if [ ! -z "$chim_cov" ]; then
        printf "$gene_id""\t""$TE_ID""\t""$chim_cov""\n" >> "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct
      fi
    fi
done <<< "$all_chim_genes"

genes=$(cut -f1 "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct | sed 's/_/\t/g' | cut -f2 | sort | uniq)
while IFS= read -r line; do
	TE_ID=$(grep "$line" "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct | cut -f2 | head -1)
	total_cov=$(grep "$line" "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct | awk '{sum+=$3;}END{print sum;}')
	isoforms=$(grep "$line" "$PROJECT"/"$sample"/tmp_dir/chimTEs_raw.ct | cut -f1 | sed -z 's/\n/,/g;s/,$/\n/')
	printf "$line""\t""$TE_ID""\t""$total_cov""\t""$isoforms""\n" >> "$PROJECT"/"$sample"/tmp_dir/"$sample"_chimTEs_final.ct
done <<< "$genes"
echo -e "${GREEN}DONE!!${NC}\n"
