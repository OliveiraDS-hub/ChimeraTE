#!/bin/bash
set -e

declare -a StringArray=("TE-initiated-UP" "TE-initiated-5UTR" "TE-exonized" "TE-terminated-3UTR" "TE-terminated-DOWN")

merging () {
raw_results=$(ls -1 "$PROJECT"/tmp/ | grep "$val" | wc -l)
if [[ "$raw_results" > 1 ]]; then
  cat "$PROJECT"/tmp/"$val"*.ct | awk '{print $1,$3}' | sed 's/ /@/g' > "$PROJECT"/tmp/"$val".lst
  freq=$(sort "$PROJECT"/tmp/"$val".lst | uniq -c | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | awk '{print $2,$1}')
  MATCH=$(awk -v number="$REP" '$2 >= number' <<< "$freq" | sed 's/_/\t/g' | awk '{print $1}')
  if [[ ! -z "$MATCH" ]]; then
    merged=$(cat "$PROJECT"/tmp/"$val"*ct)
    while IFS= read -r line;do
      TE=$(sed 's/@/\t/g' <<< "$line" | cut -f2)
      gene=$(sed 's/@/\t/g' <<< "$line" | grep "$TE" | cut -f1)
      grep "$gene" <<< "$merged" | grep "$TE" >> "$PROJECT"/tmp/"$val"_merged.lst
    done <<< "$MATCH"
  fi
fi
}

replicated () {
replicated=$(cut -f1 "$PROJECT"/tmp/"$val"_merged.lst | sort | uniq)

while IFS= read -r gene_id; do
  GENE_strd=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | cut -f2 | head -1)
  TE_IDs=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | sort -Vr -k5,5 | cut -f3 | head -1)
  TE_strd=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | sort -Vr -k5,5 | cut -f4 | head -1)
  TE_pos=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | sort -Vr -k5,5 | cut -f7 | head -1)
  EXON_pos=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | sort -Vr -k5,5 | cut -f6 | head -1)
  chim_reads=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | sort -Vr -k5,5 | cut -f5 | head -1)
  if [[ "$chim_reads" -ge "$COV" ]]; then
    printf "$gene_id""\t""$GENE_strd""\t""$TE_IDs""\t""$TE_strd""\t""$TE_pos""\t""$EXON_pos""\t""$chim_reads""\n" >> "$PROJECT"/"$val"_final.ct
  fi
done <<< "$replicated"
}

for val in ${StringArray[@]}; do
  echo -ne "Checking replicability of $val transcripts..."
  merging
  if [[ -s "$PROJECT"/tmp/"$val"_merged.lst ]]; then
    replicated
    total=$(wc -l "$PROJECT"/"$val"_final.ct | awk '{print $1}'); echo -e "${GREEN}\t==> DONE!!${NC}"; echo -e "It has been found "$total" $val transcripts in at least $REP replicates"
    echo -e "$val transcripts analysis has been finished! Check it out the result in ====> output/"$val"_final.ct\n"
    sed -i 1i"gene_id""\t""gene_strand""\t""TE_family""\t""TE_strand""\t""TE_position""\t""Exon-position""\t""chimeric_reads" "$PROJECT"/"$val"_final.ct
    rm "$PROJECT"/tmp/"$val"_merged.lst "$PROJECT"/tmp/"$val".lst; else
    echo -e "\n${YELLOW}WARNING${NC}: There is no $val transcripts found in at least "$REP" replicates!!\n"
  fi
done
