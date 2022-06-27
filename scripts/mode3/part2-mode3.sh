#!/bin/bash

set -e

skip_it () {
echo "Chimeric reads coverage..."

isoforms=$(cut -f1 "$OUTPUT"/"$sample"/all_candidates.tsv | sort | uniq)

while IFS= read -r IDs; do
        reads=$(grep "$IDs" "$OUTPUT"/"$sample"/all_candidates.tsv | cut -f2)
        while IFS= read -r id; do
                mate=$(grep '/1' <<< "$id")
                        if [ ! -z "$mate" ]; then
                                oppsit=$(sed 's|/1|/2|g' <<< "$mate")
                                chim_read=$(grep "$oppsit" "$OUTPUT"/"$sample"/transcripts_chim.lst | grep -w "$IDs")
                                else
                                mate=$(grep '/2' <<< "$id")
                                oppsit=$(sed 's|/2|/1|g' <<< "$mate")
                                chim_read=$(grep "$oppsit" "$OUTPUT"/"$sample"/transcripts_chim.lst | grep -w "$IDs")
                        fi
                        if [ ! -z "$chim_read" ]; then
                                TE_chim=$(grep -w "$id" "$OUTPUT"/"$sample"/all_candidates.tsv | cut -f4)
                                isoform=$(cut -f1 <<< "$chim_read")
                                printf "$isoform""\t""$TE_chim""\n" >> "$OUTPUT"/"$sample"/coverage.rm
                                chim_read=""
                        fi
        done <<< "$reads"

        if [ -s "$OUTPUT"/"$sample"/coverage.rm ]; then
                chim_isoform=$(cut -f1 "$OUTPUT"/"$sample"/coverage.rm | head -1)
                TE=$(cut -f2 "$OUTPUT"/"$sample"/coverage.rm | sort | uniq -c | sort -V -k 1,1 | awk '{print $2,$1}' | sed 's/ /\t/g' | head -1)
                printf "$chim_isoform""\t""$TE""\n" >> "$OUTPUT"/"$sample"/chimeric_transcripts.ct
                rm "$OUTPUT"/"$sample"/coverage.rm
        fi

done <<< "$isoforms"

echo "done"


cut -f1 "$OUTPUT"/"$sample"/chimeric_transcripts.ct > "$OUTPUT"/"$sample"/chim_trinity_IDs.lst
seqtk subseq "$OUTPUT"/"$sample"/Trinity.fasta "$OUTPUT"/"$sample"/chim_trinity_IDs.lst > "$OUTPUT"/"$sample"/chim_trinity.fa

makeblastdb -in "$TRANSCRIPTS" -dbtype nucl -out "$OUTPUT"/"$sample"/transcripts_db
blastn -query "$OUTPUT"/"$sample"/chim_trinity.fa -db "$OUTPUT"/"$sample"/transcripts_db -outfmt "6 qseqid sseqid length pident gaps mismatch qlen slen qstart qend sstart send evalue bitscore" -num_threads "$THREADS" > "$OUTPUT"/"$sample"/blast-result.tsv

isoform_IDs=$(cut -f1 "$OUTPUT"/"$sample"/blast-result.tsv | sort | uniq)

while IFS= read -r line; do
        isoform_length=$(grep "$line" "$OUTPUT"/"$sample"/blast-result.tsv | cut -f7 | head -1)
        best_hit_ID=$(grep "$line" "$OUTPUT"/"$sample"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f2 | head -1)
        best_hit_length=$(grep "$line" "$OUTPUT"/"$sample"/blast-result.tsv | grep "$best_hit_ID" | awk '{sum+=$3;} END{print sum;}')
        identity=$(grep "$line" "$OUTPUT"/"$sample"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f4 | head -1)
        ref_len=$(grep "$line" "$OUTPUT"/"$sample"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f8 | head -1)
        chim_len=$(grep "$line" "$OUTPUT"/"$sample"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f7 | head -1)
        transcript_gene=$(sed 's/_/\t/g' <<< "$best_hit_ID")
        perc=$(( best_hit_length*100/isoform_length ))
                if [ "$perc" -gt "$LENGTH" ]; then
                        printf "$line""\n" >> "$OUTPUT"/"$sample"/IDs_isoforms.lst
                        printf "$best_hit_ID""\n" >> "$OUTPUT"/"$sample"/IDs_genes.lst
                        printf "$line""\t""$transcript_gene""\t""$identity""\t""$chim_len""\t""$ref_len""\t""$best_hit_length""\n" >> "$OUTPUT"/"$sample"/blast_matches.tsv
                fi
done <<< "$isoform_IDs"

paste "$OUTPUT"/"$sample"/IDs_isoforms.lst "$OUTPUT"/"$sample"/IDs_genes.lst > "$OUTPUT"/"$sample"/chim_IDs.lst
sed -i 's/^/>/g' "$OUTPUT"/"$sample"/IDs_genes.lst

seqtk subseq "$OUTPUT"/"$sample"/chim_trinity.fa "$OUTPUT"/"$sample"/IDs_isoforms.lst > "$OUTPUT"/"$sample"/gene_homology.fa

python3 - << END
fasta= open('$OUTPUT/$sample/gene_homology.fa')
newnames= open('$OUTPUT/$sample/IDs_genes.lst')
newfasta= open('$OUTPUT/$sample/gene_homology-newIDs.fa', 'w')

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

rm "$OUTPUT"/"$sample"/IDs_isoforms.lst "$OUTPUT"/"$sample"/IDs_genes.lst
}

set +e

isoforms=$(cut -f1 "$OUTPUT"/"$sample"/chimeric_transcripts.ct)

while IFS= read -r line; do
        TE_cov=$(grep -w "$line" "$OUTPUT"/"$sample"/chimeric_transcripts.ct | cut -f2,3)
        blast_result=$(grep -w "$line" "$OUTPUT"/"$sample"/blast_matches.tsv)
        if [[ ! -z "$blast_result" ]]; then
                printf "$blast_result""\t""$TE_cov""\n" >> "$OUTPUT"/"$sample"/output.tsv
        fi
done <<< "$isoforms"


