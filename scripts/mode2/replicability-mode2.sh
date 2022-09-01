#!/bin/bash
set -e

declare -a StringArray=("output-final")

merging () {
raw_results=$(ls -1 "$PROJECT"/tmp/ | grep "$val" | wc -l)
if [[ "$raw_results" > 1 ]]; then
  cat "$PROJECT"/tmp/*"$val"*.ct | awk '{print $1,$2}' | sed 's/ /@/g' > "$PROJECT"/tmp/"$val".lst
  freq=$(sort "$PROJECT"/tmp/"$val".lst | uniq -c | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | awk '{print $2,$1}')
  MATCH=$(awk -v number="$REP" '$2 >= number' <<< "$freq" | sed 's/_/\t/g' | awk '{print $1}')
  if [[ ! -z "$MATCH" ]]; then
    merged=$(cat "$PROJECT"/tmp/*"$val"*ct)

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
  TE_IDs=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | cut -f2 | head -1)
  chim_reads=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | cut -f3)
  freq=$(wc -l <<< "$chim_reads" | awk '{print $1}')
  sum=$(awk '{sum+=$1;} END{print sum;}' <<< "$chim_reads")
  avg_reads=$(awk -v var1=$sum -v var2=$freq 'BEGIN { print  ( var1 / var2 ) }' | awk '{printf("%.4f\n", $1)}')

  isoforms=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | cut -f4 | sed -z 's/\n/,/g;s/,$/\n/;s/,/\n/g' | sort | uniq -c | awk -v number="$REP" '$1 >= number' | awk '{print $2}' | sed -z 's/\n/,/g;s/,$/\n/')
  expression=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | cut -f5)
  freq_exp=$(wc -l <<< "$expression" | awk '{print $1}')
  sum_exp=$(awk '{sum+=$1;} END{print sum;}' <<< "$expression")
  avg_fpkm=$(awk -v var1=$sum_exp -v var2=$freq_exp 'BEGIN { print  ( var1 / var2 ) }' | awk '{printf("%.4f\n", $1)}')

  if (( $(echo "$avg_fpkm $COV" | awk '{print ($1 >= $2)}') )); then
	    printf "$gene_id""\t""$TE_IDs""\t""$avg_reads""\t""$isoforms""\t""$avg_fpkm""\n" >> "$PROJECT"/chimTE-final-chimreads.ct
	fi

done <<< "$replicated"
}

for val in ${StringArray[@]}; do
  merging
  if [[ -s "$PROJECT"/tmp/"$val"_merged.lst ]]; then
    replicated
    if [[ -f "$PROJECT"/chimTE-final-chimreads.ct && -z "$ASSEMLY" ]]; then
      total=$(wc -l "$PROJECT"/chimTE-final-chimreads.ct | awk '{print $1}'); echo -e "It has been found "$total" chimeric transcripts only with chimeric reads evidence"
      echo -e "The analysis has been finished! Check it out the result in ====> "$PROJECT"/chimTE-final-chimreads.ct\n"; else
      echo -e "The analysis did not find any chimeric transcripts based on chimeric reads evidence!"
    fi
    rm "$PROJECT"/tmp/*"$val"*.lst
  fi
done
