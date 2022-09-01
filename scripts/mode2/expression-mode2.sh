#!/bin/bash

echo -e "Calculating chimeric reads coverage..."

gene_list=$(cut -f1 "$PROJECT"/"$sample"/tmp_dir/"$sample"_chimTEs_final.ct)
transcripts_list=$(cut -f4 "$PROJECT"/"$sample"/tmp_dir/"$sample"_chimTEs_final.ct | sed 's/,/\n/g')
cut -f4 "$PROJECT"/"$sample"/tmp_dir/"$sample"_chimTEs_final.ct | sed 's/,/\n/g' > "$PROJECT"/"$sample"/tmp_dir/transcripts_list.lst

while IFS= read -r gene_id; do
  spc_transcripts=$(grep -F "$gene_id" "$PROJECT"/"$sample"/tmp_dir/transcripts_list.lst)

  while IFS= read -r transcript_id; do
    spc_fpkm=$(grep "$transcript_id" "$ALN"/fpkm_counts/expressed_transcripts.xprs | awk '{print $13}' | awk '{printf("%.4f\n", $1)}')
    all_fpkm+="${spc_fpkm}\n"
    all_Rnumber=$(echo -e "$all_fpkm")
  done <<< "$spc_transcripts"

  freq=$(wc -l <<< "$all_Rnumber" | awk '{print $1}')
  sum=$(awk '{sum+=$1;} END{print sum;}' <<< "$all_Rnumber")
  avg_fpkm=$(awk -v var1=$sum -v var2=$freq 'BEGIN { print  ( var1 / var2 ) }' | awk '{printf("%.4f\n", $1)}')

  if (( $(echo "$avg_fpkm $FPKM" | awk '{print ($1 >= $2)}') )); then
    result=$(grep -w "$gene_id" "$PROJECT"/"$sample"/tmp_dir/"$sample"_chimTEs_final.ct)
    paste <(echo "$result") <(echo "$avg_fpkm") --delimiters '\t' >> "$PROJECT"/"$sample"/"$sample"_output-final.ct
  fi

  all_fpkm=""

done <<< "$gene_list"
echo -e "${GREEN}DONE!!${NC}\n"
