#!/bin/bash

set +e
echo -ne "Chimeric reads coverage..."

isoforms=$(cut -f1 "$TRINITY_OUT"/all_candidates.tsv | sort | uniq)

while IFS= read -r IDs; do
        reads=$(grep "$IDs" "$TRINITY_OUT"/all_candidates.tsv | cut -f2)
        while IFS= read -r id; do
                mate=$(grep '/1' <<< "$id")
                        if [ ! -z "$mate" ]; then
                                oppsit=$(sed 's|/1|/2|g' <<< "$mate")
                                chim_read=$(grep "$oppsit" "$TRINITY_OUT"/transcripts_chim.lst | grep -w "$IDs")
                                else
                                mate=$(grep '/2' <<< "$id")
                                oppsit=$(sed 's|/2|/1|g' <<< "$mate")
                                chim_read=$(grep "$oppsit" "$TRINITY_OUT"/transcripts_chim.lst | grep -w "$IDs")
                        fi
                        if [ ! -z "$chim_read" ]; then
                                TE_chim=$(grep -w "$id" "$TRINITY_OUT"/all_candidates.tsv | cut -f4)
                                isoform=$(cut -f1 <<< "$chim_read")
                                printf "$isoform""\t""$TE_chim""\n" >> "$TRINITY_OUT"/coverage.rm
                                chim_read=""
                        fi
        done <<< "$reads"

        if [ -s "$TRINITY_OUT"/coverage.rm ]; then
                chim_isoform=$(cut -f1 "$TRINITY_OUT"/coverage.rm | head -1)
                TE=$(cut -f2 "$TRINITY_OUT"/coverage.rm | sort | uniq -c | sort -V -k 1,1 | awk '{print $2,$1}' | sed 's/ /\t/g' | head -1)
                printf "$chim_isoform""\t""$TE""\n" >> "$TRINITY_OUT"/chimeric_transcripts.ct
                rm "$TRINITY_OUT"/coverage.rm
        fi

done <<< "$isoforms"
echo -e "\t${GREEN}===> DONE!${NC}\n"

set -e
echo -ne "Performing blast to identify transcripts..."
cut -f1 "$TRINITY_OUT"/chimeric_transcripts.ct > "$TRINITY_OUT"/chim_trinity_IDs.lst
seqtk subseq "$TRINITY_OUT"/Trinity.fasta "$TRINITY_OUT"/chim_trinity_IDs.lst > "$TRINITY_OUT"/chim_trinity.fa

makeblastdb -in "$TRANSCRIPTS_FASTA" -dbtype nucl -out "$TRINITY_OUT"/transcripts_db >> /dev/null 2>&1
blastn -query "$TRINITY_OUT"/chim_trinity.fa -db "$TRINITY_OUT"/transcripts_db -outfmt "6 qseqid sseqid length pident gaps mismatch qlen slen qstart qend sstart send evalue bitscore" -num_threads "$THREADS" > "$TRINITY_OUT"/blast-result.tsv
echo -e "\t${GREEN}===> DONE!${NC}\n"

echo -ne "Recovering the best matches..."
isoform_IDs=$(cut -f1 "$TRINITY_OUT"/blast-result.tsv | sort | uniq)

while IFS= read -r line; do
        isoform_length=$(grep "$line" "$TRINITY_OUT"/blast-result.tsv | cut -f7 | head -1)
        best_hit_ID=$(grep "$line" "$TRINITY_OUT"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f2 | head -1)
        best_hit_length=$(grep "$line" "$TRINITY_OUT"/blast-result.tsv | grep "$best_hit_ID" | awk '{sum+=$3;} END{print sum;}')
        identity=$(grep "$line" "$TRINITY_OUT"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f4 | head -1)
        ref_len=$(grep "$line" "$TRINITY_OUT"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f8 | head -1)
        chim_len=$(grep "$line" "$TRINITY_OUT"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f7 | head -1)
        transcript_gene=$(sed 's/_/\t/g' <<< "$best_hit_ID")
        perc=$(( best_hit_length*100/isoform_length ))
                if [ "$perc" -gt "$LENGTH" ]; then
                        printf "$line""\n" >> "$TRINITY_OUT"/IDs_isoforms.lst
                        printf "$best_hit_ID""\n" >> "$TRINITY_OUT"/IDs_genes.lst
                        printf "$line""\t""$transcript_gene""\t""$identity""\t""$chim_len""\t""$ref_len""\t""$best_hit_length""\n" >> "$TRINITY_OUT"/blast_matches.tsv
                fi
done <<< "$isoform_IDs"
echo -e "\t${GREEN}===> DONE!${NC}\n"

echo -ne "Creating output file..."
paste "$TRINITY_OUT"/IDs_isoforms.lst "$TRINITY_OUT"/IDs_genes.lst > "$TRINITY_OUT"/chim_IDs.lst
sed -i 's/^/>/g' "$TRINITY_OUT"/IDs_genes.lst

seqtk subseq "$TRINITY_OUT"/chim_trinity.fa "$TRINITY_OUT"/IDs_isoforms.lst > "$TRINITY_OUT"/gene_homology.fa

python3 - << END
fasta= open('$TRINITY_OUT/gene_homology.fa')
newnames= open('$TRINITY_OUT/IDs_genes.lst')
newfasta= open('$TRINITY_OUT/gene_homology-newIDs.fa', 'w')

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

rm "$TRINITY_OUT"/IDs_isoforms.lst "$TRINITY_OUT"/IDs_genes.lst

isoforms=$(cut -f1 "$TRINITY_OUT"/chimeric_transcripts.ct)
set +e
while IFS= read -r line; do
        TE_cov=$(grep -w "$line" "$TRINITY_OUT"/chimeric_transcripts.ct | cut -f2,3)
        blast_result=$(grep -w "$line" "$TRINITY_OUT"/blast_matches.tsv)
        if [[ ! -z "$blast_result" ]]; then
                printf "$blast_result""\t""$TE_cov""\n" >> "$TRINITY_OUT"/"$sample"_output.tsv
        fi
done <<< "$isoforms"
echo -e "\t${GREEN}===> DONE!${NC}\n  "

#echo -e "Quantifying transcripts expression..."
#mkdir  "$TRINITY_OUT"/expression 2>/dev/null
#align_and_estimate_abundance.pl --transcripts "$TRINITY_OUT"/Trinity.fasta --seqType fq --left "$mate1_rep" --right "$mate2_rep" --est_method RSEM --aln_method bowtie2 --SS_lib_type "$STRANDED" --thread_count "$THREADS" --gene_trans_map  "$TRINITY_OUT"/Trinity.fasta.gene_trans_map --output_dir "$TRINITY_OUT"/expression --prep_reference
#echo -ne "\tDONE!"
