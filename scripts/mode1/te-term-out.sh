#!/bin/bash

# Merging UTR cov and Upstream cov

merging_3UTR () {
cut -f1 "$PROJECT"/TE-terminated.ct | sort | uniq > "$PROJECT"/TE-terminated_genes.lst

if [ -a "$PROJECT"/TE-terminated_UTRs.ct ]; then
	if [ -a "$PROJECT"/TE-terminated.ct ]; then
		while IFS= read -r line
		do
			UP_info=$(grep -w "$line" "$PROJECT"/TE-terminated.ct)
			UP_gene_ID=$(cut -f1 <<< "$UP_info" | sort | uniq)
			UP_TE_number=$(grep -F -w "$line" <<< "$UP_info" | cut -f3 | sort | uniq | wc -l)
			UP_cov=$(grep -F -w "$line" <<< "$UP_info" | cut -f5)
				if [ "$UP_TE_number" -gt 1 ]; then
					UP_cov=$(awk '{ total += $1 } END { print total/NR }' <<< "$UP_cov")
				fi
			UP_TE_IDs=$(grep -F -w "$line" <<< "$UP_info" | cut -f3 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
			UP_TE_pos=$(grep -F -w "$line" <<< "$UP_info" | cut -f7 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')

			UTR_info=$(grep -w "$line" "$PROJECT"/TE-terminated_UTRs.ct)
				if [ ! -z "$UTR_info" ]; then
					UTR_TE_number=$(grep -F -w "$line" <<< "$UTR_info" | cut -f3 | sort | uniq | wc -l)
					UTR_cov=$(grep -F -w "$line" <<< "$UTR_info" | cut -f5)
					UTR_TE_IDs=$(grep -F -w "$line" <<< "$UTR_info" | cut -f3 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
					UTR_number=$(grep -F -w "$line" <<< "$UTR_info" | cut -f6 | sort | uniq | wc -l)
					nTEs=$(grep -F -w "$line" <<< "$UTR_info" | cut -f3 | sort | uniq | wc -l)
						if [ "$nTEs" -gt 1 ]; then
							UTR_cov=$(awk '{ total += $1 } END { print total/NR }' <<< "$UTR_cov"); else
							UTR_cov=$(head -1 <<< "$UTR_cov")
						fi
					UTR_pos=$(grep -F -w "$line" <<< "$UTR_info" | cut -f6 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
					UTR_TE_pos=$(grep -F -w "$line" <<< "$UTR_info" | cut -f7 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g'); else

					UTR_TE_number="0"
					UTR_cov="0"
					UTR_TE_IDs="-"
					UTR_number="0"
					UTR_pos="-"
					UTR_TE_pos="-"
				fi

				printf "$UP_gene_ID""\t""$UP_TE_number""\t""$UP_TE_IDs""\t""$UP_TE_pos""\t""$UP_cov""\t""$UTR_TE_IDs""\t""$UTR_TE_pos""\t""$UTR_number""\t""$UTR_pos""\t""$UTR_cov""\n" >> "$PROJECT"/TE-terminated_merging.ct

		done < "$PROJECT"/TE-terminated_genes.lst
	else
	echo "te-terminated.ct doesnt exist"
	fi
else
echo "5utr doesnt exist"
fi

IDs=$(cut -f1 "$PROJECT"/TE-terminated.ct | sort | uniq)
o_UTR=$(grep -v -F -w "$IDs" "$PROJECT"/TE-terminated_UTRs.ct | cut -f1 | sort | uniq)

while read -r line
do
	UP_TE_number="0"
	UP_TE_IDs="-"
	UP_TE_pos="-"
	UP_cov="0"

	UTR_TE_IDs=$(grep -F -w "$line" "$PROJECT"/TE-terminated_UTRs.ct | cut -f3 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
	UTR_TE_pos=$(grep -F -w "$line" "$PROJECT"/TE-terminated_UTRs.ct | cut -f7 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
	UTR_number=$(grep -F -w "$line" "$PROJECT"/TE-terminated_UTRs.ct | cut -f6 | sort | uniq | wc -l)
	UTR_pos=$(grep -F -w "$line" "$PROJECT"/TE-terminated_UTRs.ct | cut -f6 | sort | uniq | tr '\n' ' ' | sed 's/ /; /g' | sed 's/..$//g')
	UTR_cov=$(grep -F -w "$line" "$PROJECT"/TE-terminated_UTRs.ct | cut -f5)
	UTR_TE_number=$(grep -F -w "$line" "$PROJECT"/TE-terminated_UTRs.ct | cut -f3 | sort | uniq | wc -l)
		if [ "$UTR_number" -gt 1 ]; then
			UTR_cov=$(awk '{ total += $1 } END { print total/NR }' <<< "$UTR_cov")
		fi
	printf "$line""\t""$UP_TE_number""\t""$UP_TE_IDs""\t""$UP_TE_pos""\t""$UP_cov""\t""$UTR_TE_IDs""\t""$UTR_TE_pos""\t""$UTR_number""\t""$UTR_pos""\t""$UTR_cov""\n" >> "$PROJECT"/TE-terminated_merging.ct
done <<< "$o_UTR"

cutoff="2"
while IFS= read -r line
do
	cov=$(cut -f5 <<< "$line")
	cov_utr=$(cut -f10 <<< "$line")
	if (( $(bc <<< "$cov >= $cutoff" ) || $(bc <<< "$cov_utr >= $cutoff" ) )); then
		echo "$line" >> "$PROJECT"/TE-terminated_final.ct
	fi
done < "$PROJECT"/TE-terminated_merging.ct

sed -i '1i\Gene_ID\tTEs_upstream\tTE_IDs_upstream\tTE_upstrem_position\tUpstream_cov\tTE_IDs_5UTR\tTE_position_5UTR\t5UTR_number\t5UTR_position\t5UTR_cov' "$PROJECT"/TE-terminated_final.ct
}
if [ ! -z $UTR ] ; then
        merging_3UTR
fi
