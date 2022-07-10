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
  freq=$(wc -l <<< "$chim_reads")
  sum=$(perl -lne '{chomp;$x+=$_}END{print $x}' <<< $chim_reads)
  avg_chim_reads=$(echo "scale=4; $sum/$freq" | bc)
  avg_2f=$(LC_NUMERIC="en_US.UTF-8" printf '%0.2f\n' "$avg_chim_reads")

  isoforms=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | cut -f4 | sed -z 's/\n/,/g;s/,$/\n/;s/,/\n/g' | sort | uniq -c | awk -v number="$REP" '$1 >= number' | awk '{print $2}' | sed -z 's/\n/,/g;s/,$/\n/')

  expression=$(grep -F -w "$gene_id" "$PROJECT"/tmp/"$val"_merged.lst | cut -f5)
  freq_exp=$(wc -l <<< "$expression")
  sum_exp=$(perl -lne '{chomp;$x+=$_}END{print $x}' <<< $expression)

  avg_fpkm=$(echo "scale=4; $sum_exp/$freq_exp" | bc)
  avg_fpkm_2f=$(LC_NUMERIC="en_US.UTF-8" printf '%0.2f\n' "$avg_fpkm")

  if [ 1 -eq "$(echo "${avg_2f} >= ${COV}" | bc)" ]; then
    printf "$gene_id""\t""$TE_IDs""\t""$avg_2f""\t""$isoforms""\t""$avg_fpkm_2f""\n" >> "$PROJECT"/chimeraTE_final.ct
  fi
done <<< "$replicated"
}

for val in ${StringArray[@]}; do
  echo -e "############### Checking replicability of chimeric transcripts ###############\n"
  merging
  if [[ -s "$PROJECT"/tmp/"$val"_merged.lst ]]; then
    replicated
    total=$(wc -l "$PROJECT"/chimeraTE_final.ct | awk '{print $1}'); echo -e "It has been found "$total" chimeric transcripts in at least $REP replicates!"
    echo -e "$val transcripts analysis has been finished! Check it out the result in ====> "$val"_final.ct\n"
    rm "$PROJECT"/tmp/*"$val"*.lst
  else
    echo -e "There is no $val transcripts found in at least "$REP" replicates!!\n"
  fi
done
