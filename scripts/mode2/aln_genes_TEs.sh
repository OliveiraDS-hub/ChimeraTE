#!/bin/bash

max=$(grep -c '>' "$TE_FASTA"); i=0
grep '>' "$TE_FASTA" > "$PROJECT"/"$sample"/tmp_dir/TE_IDs.lst

while [ $i -ne $max ]; do
  while IFS= read -r line; do
    i=$(($i+1))
    echo "$line&$i" >> "$PROJECT"/"$sample"/tmp_dir/new_IDs.lst
  done < "$PROJECT"/"$sample"/tmp_dir/TE_IDs.lst
done

python3 - << END
fasta= open('$TE_FASTA')
newnames= open('$PROJECT/$sample/tmp_dir/new_IDs.lst')
newfasta= open('$PROJECT/$sample/tmp_dir/TE_seq_mod.fa', 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()
END

TE_file_mod="$PROJECT/$sample/tmp_dir/TE_seq_mod.fa"
TE_file=$(basename "$TE_FASTA" | sed 's/.fasta//g' | sed 's/.fa//g')

if [ -f "$PROJECT"/indexes/"$TE_file".1.bt2 ]; then
echo -e "\nBowtie2 index for $TE_file has been found"; else
echo -e "Creating Bowtie2 index with $TE_file..."
bowtie2-build "$TE_file_mod" "$PROJECT"/indexes/"$TE_file" --threads "$THREADS" >> /dev/null 2>&1; echo -e "${GREEN}DONE!!${NC}\n"
fi

echo -e "Running TEs alignment..."
bowtie2 -x "$PROJECT"/indexes/"$TE_file" -1 "$mate1" -2 "$mate2" -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p "$THREADS" -S "$ALN"/tes.sam
echo -e "${GREEN}DONE!!${NC}\n"

transcripts_file=$(basename "$TRANSCRIPTS_FASTA" | sed 's/.fasta//g' | sed 's/.fa//g')
if [ -f "$PROJECT"/indexes/"$transcripts_file".1.bt2 ]; then
echo -e "\nBowtie2 index for $PROJECT has been found"; else
echo -e "Creating Bowtie2 index with $TE_file..."
bowtie2-build "$TRANSCRIPTS_FASTA" "$PROJECT"/indexes/"$transcripts_file" --threads "$THREADS" >> /dev/null 2>&1; echo -e "${GREEN}DONE!!${NC}\n"
fi

echo -e "Running genes alignment..."
bowtie2 -x "$PROJECT"/indexes/"$transcripts_file" -1 "$mate1" -2 "$mate2" -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p "$THREADS" -S "$ALN"/genes.sam
echo -e "${GREEN}DONE!!${NC}\n"

echo -e "Converting sam to bed..."
for SAMF in "$ALN"/*.sam; do
	output=$(sed 's|.*/||g' <<< "$SAMF" | sed 's/.sam//g')
	samtools view -@ "$THREADS" -bS "$SAMF" > "$ALN"/"$output".bam 2>/dev/null; rm "$SAMF"
	bedtools bamtobed -i  "$ALN"/"$output".bam > "$ALN"/"$output"_total.bed
done
awk '{gsub("&.*","",$1)}1' OFS='\t' "$ALN"/tes_total.bed > "$ALN"/tmp_file.tmp && mv "$ALN"/tmp_file.tmp "$ALN"/tes_total.bed
echo -e "${GREEN}DONE!!${NC}\n"

echo -e "Calculating FPKM to the transcripts..."
express -o "$ALN"/fpkm_counts -O 1 --output-align-prob --no-bias-correct --"$STRANDED" "$TRANSCRIPTS_FASTA" "$ALN"/genes.bam >> /dev/null 2>&1; rm "$ALN"/fpkm_counts/hits.1.prob.bam >> /dev/null 2>&1
echo -e "${GREEN}DONE!!${NC}\n"

echo -e "Removing transcripts with FPKM lower than 1..."
awk '$13 >= "1"' "$ALN"/fpkm_counts/results.xprs > "$ALN"/fpkm_counts/expressed_transcripts.xprs
awk '$13 < "1"' "$ALN"/fpkm_counts/results.xprs | cut -f2 | sed '1d' > "$ALN"/non_exp_transcripts.lst
LC_ALL=C fgrep -w -v -f "$ALN"/non_exp_transcripts.lst "$ALN"/genes_total.bed > "$ALN"/genes_total_expressed.bed
echo -e "${GREEN}DONE!!${NC}\n"

echo -e "Identifying chimeric reads..."

cut -f4 "$ALN"/tes_total.bed | grep '/1' | sed 's|/1||' > "$ALN"/TEs_R1_IDs.lst
cut -f4 "$ALN"/tes_total.bed | grep '/2' | sed 's|/2||' > "$ALN"/TEs_R2_IDs.lst
cut -f4 "$ALN"/genes_total_expressed.bed | grep '/1' | sed 's|/1||' > "$ALN"/genes_R1_IDs.lst
cut -f4 "$ALN"/genes_total_expressed.bed | grep '/2' | sed 's|/2||' > "$ALN"/genes_R2_IDs.lst

LC_ALL=C fgrep -w -f "$ALN"/TEs_R1_IDs.lst "$ALN"/genes_R2_IDs.lst | sort | uniq > "$READS"/all_chim_IDs_1.lst
LC_ALL=C fgrep -w -f "$ALN"/TEs_R2_IDs.lst "$ALN"/genes_R1_IDs.lst | sort | uniq > "$READS"/all_chim_IDs_2.lst

if [[ -f "$READS"/all_chim_IDs_1.lst || -f "$READS"/all_chim_IDs_2.lst ]]; then
	cat "$READS"/all_chim_IDs_1.lst "$READS"/all_chim_IDs_2.lst | sort | uniq > "$READS"/all_chim_IDs.lst; else
	echo "There is no chimeric read pairs! The analysis cannot cotinue. Exiting..." && exit 1
fi

LC_ALL=C fgrep -w -f "$READS"/all_chim_IDs.lst "$ALN"/tes_total.bed | sed 's|/2||; s|/1||' > "$READS"/TEs_chim.bed
LC_ALL=C fgrep -w -f "$READS"/all_chim_IDs.lst "$ALN"/genes_total_expressed.bed | sed 's|/2||; s|/1||' > "$READS"/genes_chim.bed
echo -e "${GREEN}DONE!!${NC}\n"
