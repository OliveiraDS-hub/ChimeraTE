#!/bin/bash

cat "$PROJECT"/tmp/*.ct | cut -f1,2 | sed 's/\t/_/g' > "$PROJECT"/tmp/chimeras.lst
#replicates=$(cut -f1 $1 | sed 's|.*/||g; s|/.*/||g; s/.fastq//g; s/.fq//g; s/.gz//g; s/_R1//g')

#while IFS= read -r ID; do
	#cut -f1,2 projects/"$ID"/"$ID"_chimTEs_final.ct | sed 's/\t/_/g' > projects/replicates/"$ID".lst
	#cp projects/"$ID"/"$ID"_chimTEs_final.ct projects/replicates/
#done <<< "$replicates"

#cat projects/replicates/*final.ct > projects/replicates/all_chimTE_results.ct
#cat projects/replicates/*.lst > projects/replicates/all_IDs.lst

#total_rep=$(wc -l <<< "$replicates")
#echo "$total_rep"

freq=$(sort "$PROJECT"/tmp/chimeras.lst | uniq -c | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | awk '{print $2,$1}')

wc -l <<< "$freq"
#MATCH=$(awk -v number="$total_rep" '$2 == number' <<< "$freq" | sed 's/_/\t/g' | awk '{print $1}')

pause () {
head <<< "$MATCH"
echo "Merging chimeras from all replicates"
while IFS= read -r line; do
	GENE_ID=$(grep -F -w "$line" projects/replicates/all_chimTE_results.ct | cut -f1 | head -1)
	TE_IDs=$(grep -F -w "$line" projects/replicates/all_chimTE_results.ct | cut -f2 | head -1)
	cov=$(grep -F -w "$line" projects/replicates/all_chimTE_results.ct | cut -f3)
	avg=$(awk '{ total += $1 } END { print total/NR }' <<< "$cov")
	#if [ "$avg" -ge "$cutoff" ]; then
	if (( $(echo "$avg >= $cutoff" |bc -l) )); then
        	printf "$GENE_ID""\t""$TE_IDs""\t""$avg""\n" >> projects/replicates/$2_chimTEs_merging.ct
	fi
done <<< "$MATCH"
}
