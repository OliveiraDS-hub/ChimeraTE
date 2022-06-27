#!/bin/bash

tudo () {
echo -e "Calculating chimeric reads coverage...\n"

expression=$(cat "$PROJECT"/"$sample"/tmp_dir/"$sample"_chimTEs_final.ct)

gene_list=$(sed '1d' <<< "$expression" | cut -f1)
transcripts_list=$(sed '1d' <<< "$expression" | cut -f4 | sed 's/,/\n/g')
sed '1d' <<< "$expression" | cut -f4 | sed 's/,/\n/g' | sed 's/_/\t/g' > "$PROJECT"/"$sample"/tmp_dir/transcripts_list.lst

express=$(cat "$ALN"/fpkm_counts/expressed_transcripts.xprs)
min_fpkm="1"

while IFS= read -r gene_id; do
  spc_transcripts=$(grep -w "$gene_id" "$PROJECT"/"$sample"/tmp_dir/transcripts_list.lst | cut -f1)

    while IFS= read -r transcript_id; do
      spc_fpkm=$(grep "$transcript_id" "$ALN"/fpkm_counts/expressed_transcripts.xprs | awk '{print $13}')
      r_number=$(LC_NUMERIC="en_US.UTF-8" printf '%.4f' "$spc_fpkm")
      all_fpkm+="${r_number}\n"
      all_Rnumber=$(echo -e "$all_fpkm")
    done <<< "$spc_transcripts"

  freq=$(wc -l <<< "$all_Rnumber" | awk '{print $1}')
  sum=$(perl -lne '{chomp;$x+=$_}END{print $x}' <<< $all_Rnumber)
  avg_fpkm=$(echo "$sum / $freq" | bc -l)
  avg_4f=$(LC_NUMERIC="en_US.UTF-8" printf '%0.4f\n' "$avg_fpkm")
  #echo -e "$gene_id \t $avg_4f"
  all_fpkm=""

  result=$(grep -w "$gene_id" "$PROJECT"/"$sample"/tmp_dir/"$sample"_chimTEs_final.ct)
  paste <(echo "$result") <(echo "$avg_4f") --delimiters '\t' >> "$PROJECT"/"$sample"/"$sample"_output-final.ct
done <<< "$gene_list"
}

tudo

cancel () {

freq=$(wc -l sci-number.lst | awk '{print $1}')
echo "freq = $freq"

while IFS= read -r gene_id; do
  real_number=$(LC_NUMERIC="en_US.UTF-8" printf '%.4f\n' "$gene_id")
  all_Rnumber+="${real_number}\n"

done < sci-number.lst

all_Rnumber=$(echo -e "$all_Rnumber")

#echo -e "all:\n$all_Rnumber"
#echo -e "$all_Rnumber" > file.txt
#sum=$(awk '{a+=$0}END{print a}' file.txt)
sum=$(perl -lne '{chomp;$x+=$_}END{print $x}' <<< $all_Rnumber)
echo "sum = $sum"

avg_fpkm=$(echo "$sum / $freq" | bc -l)
echo "avg = $avg_fpkm"
LC_NUMERIC="en_US.UTF-8" printf '%0.4f\n' "$avg_fpkm"
}
#0.00000000004253




#for file in *
#do
#    lineInfo=`wc -l $file`
#    total="$total$lineInfo, "  # or total+="$lineInfo, "
#done
#awk '{ sum += $2 } END { if (NR > 0) print sum / NR }'




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
