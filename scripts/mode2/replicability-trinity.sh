#!/bin/bash
set -e
raw_results=$(ls -1 "$PROJECT"/tmp | grep "output.tsv")

simplyfing () {
n_files=$(ls -1 "$PROJECT"/tmp | grep "output.tsv" | wc -l)
if [[ "$n_files" > 1 ]]; then
  while IFS= read -r raw_file; do
    genes=$(cut -f3 "$PROJECT"/tmp/"$raw_file" | sort | uniq)
    while IFS= read -r ID; do
      genelist=$(grep -w "$ID" "$PROJECT"/tmp/"$raw_file")
      top=$(sort -Vr -k9,9 <<< "$genelist" | head -1)
      echo -e "$top" >> "$PROJECT"/tmp/"$raw_file"_final.lst
    done <<< "$genes"
  done <<< "$raw_results"
fi
}

replicated () {
freq=$(cat "$PROJECT"/tmp/*output.tsv_final.lst | awk '{print $3,$8}' | sed 's/ /@/g' | sort | uniq -c)
MATCH=$(awk -v number="$REP" '$1 >= number' <<< "$freq" | awk '{print $2}' | sed 's/@/\t/g' | awk '{print $1}')
raw_data=$(ls -1 "$PROJECT"/tmp/ | grep "_final.lst")

while IFS= read -r file; do
  while IFS= read -r gene_id; do
    grep -w "$gene_id" "$PROJECT"/tmp/"$file" >> "$PROJECT"/tmp/replicated.tsv
  done <<< "$MATCH"

done <<< "$raw_data"
}

final () {

replicated_genes=$(cut -f3 "$PROJECT"/tmp/replicated.tsv | sort | uniq)

while IFS= read -r gene_id; do
  chimeras=$(grep -w "$gene_id" "$PROJECT"/tmp/replicated.tsv | sort -Vr -k5,5 | head -1 | cut -f1-8)
  chim_reads=$(grep -w "$gene_id" "$PROJECT"/tmp/replicated.tsv | cut -f9)
  freq=$(grep -w "$gene_id" "$PROJECT"/tmp/replicated.tsv | wc -l | awk '{print $1}')
  sum=$(awk '{sum+=$1;} END{print sum;}' <<< "$chim_reads")
  avg_reads=$(awk -v var1=$sum -v var2=$freq 'BEGIN { print  ( var1 / var2 ) }' | awk '{printf("%.4f\n", $1)}')
  printf "$chimeras""\t""$avg_reads""\n" >> "$PROJECT"/tmp/replicated_final.tsv
done <<< "$replicated_genes"

}

double-evidence () {
  cut -f3 "$PROJECT"/tmp/replicated_final.tsv > "$PROJECT"/tmp/geneIDs-assembly.lst
  cut -f1 "$PROJECT"/chimTE-final-chimreads.ct > "$PROJECT"/tmp/geneIDs-chimTE.lst
  cat "$PROJECT"/tmp/geneIDs-* | sort | uniq -c | awk '$1 >= 2' | awk '{print $2}' > "$PROJECT"/tmp/double-evidence.lst

  while IFS= read -r gene_id; do
    chimTE=$(grep -w "$gene_id" "$PROJECT"/chimTE-final-chimreads.ct)
    isoform_ID=$(grep -w "$gene_id" "$PROJECT"/tmp/replicated_final.tsv | cut -f1)
    rest=$(grep -w "$gene_id" "$PROJECT"/tmp/replicated_final.tsv | cut -f4-9)
    echo -e "$chimTE\t$isoform_ID\t$rest" >> "$PROJECT"/chimTE-final-double-evidence.ct
  done < "$PROJECT"/tmp/double-evidence.lst

  solo_assembly=$(cat "$PROJECT"/tmp/double-evidence.lst "$PROJECT"/tmp/geneIDs-assembly.lst | sort | uniq -c | awk '$1 == "1"' | awk '{print $2}')
  while IFS= read -r gene_id; do
    grep -w "$gene_id" "$PROJECT"/tmp/replicated_final.tsv >> "$PROJECT"/chimTE-final-transcriptome.ct
  done <<< "$solo_assembly"

  solo_chimTE=$(cat "$PROJECT"/tmp/double-evidence.lst "$PROJECT"/tmp/geneIDs-chimTE.lst | sort | uniq -c | awk '$1 == "1"' | awk '{print $2}')
  while IFS= read -r gene_id; do
    grep -w "$gene_id" "$PROJECT"/chimTE-final-chimreads.ct >> "$PROJECT"/chimTE-final-chimreads-without-assembly.ct
  done <<< "$solo_chimTE"

  rm "$PROJECT"/tmp/geneIDs* "$PROJECT"/tmp/double-evidence.lst "$PROJECT"/tmp/replicated.tsv
}

simplyfing
replicated
final
double-evidence
mv "$PROJECT"/chimTE-final-chimreads.ct "$PROJECT"/tmp/

if [[ -z "$PROJECT"/chimTE-final-transcriptome.ct ]]; then
  total=$(wc -l "$PROJECT"/chimTE-final-transcriptome.ct | awk '{print $1}'); echo -e "It has been found "$total" chimeric transcripts only with transcriptome assembly"
  echo -e "Check it out the result in ====> $PROJECT/chimTE-final-transcriptome.ct"; else
  echo -e "The analysis did not find any chimeric transcripts based on transcriptome assembly!"
fi

if [[ -z "$PROJECT"/chimTE-final-chimreads-without-assembly.ct ]]; then
  total=$(wc -l "$PROJECT"/chimTE-final-chimreads-without-assembly.ct | awk '{print $1}'); echo -e "It has been found "$total" chimeric transcripts only with chimeric reads alignment"
  echo -e "Check it out the result in ====> $PROJECT/chimTE-final-chimreads-without-assembly.ct"; else
  echo -e "The analysis did not find any chimeric transcripts based on transcriptome assembly!"
fi

if [[ -z "$PROJECT"/chimTE-final-double-evidence.ct ]]; then
  total=$(wc -l "$PROJECT"/chimTE-final-double-evidence.ct | awk '{print $1}'); echo -e "It has been found "$total" chimeric transcripts based on both chimeric reads alignemnt and transcriptome assembly"
  echo -e "Check it out the result in ====> $PROJECT/chimTE-final-double-evidence.ct"; else
  echo -e "The analysis did not find any chimeric transcripts based on both chimeric reads alignemnt and transcriptome assembly!"
fi
