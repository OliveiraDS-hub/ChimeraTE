#!/bin/bash

set -a

function usage() {
   cat << help

######################################
######## ChimeraTE v1.0 mode3 ########
######################################

USAGE:

-One-replicate:

ChimTE-transcripts.sh [--mate1 <mate1.fastq.gz>] [--mate2 <mate2.fastq.gz>] [--stranded <RF or FR>] [--transcripts <transcripts.fa>] [--ref_TEs <taxonomy or library_TEs.fa>] [options]

-Multi-replicates:

ChimTE-transcripts.sh [--mate1 <mate1_replicate1.fastq.gz,mate1_replicate2.fastq.gz>] [--mate2 <mate2_replicate1.fastq.gz,mate2_replicate2.fastq.gz>] [--stranded <RF or FR>] [--transcripts <transcripts.fa>] [--ref_TEs <taxonomy or library_TEs.fa>] [options]


#Mandatory arguments:

    --mate1			mate 1 from paired-end reads

    --mate2			mate 2 from paired-end reads

    --stranded		  	indicates the order of the stranded RNA-seq (RF - reverse/forward; FR - forward/reverse)

    --transcripts		fasta file with reference transcripts

    --ref_TEs			"species" database used by RepeatMasker (flies, human, mouse, arabidopsis);
				or a built TE library in fasta format

#Optional arguments:

    --output		  	Output directory

    --threads               	Number of threads, (default: 6)

    --ram                   	RAM memory (default: 8 Gbytes)

    --TE_length                 Minimum TE length to keep it from RepeatMasker output (default: 80bp)

    --min_length          	Minimum length with homology between de novo assembled transcripts and reference transcripts (default: 80%)

    --overlap              	Minmum overlap (0.1 to 1) between read length and TE insertion (default: 0.5)

    --min_TPM             	Minimum TPM expression (default: 1)

    -h, --help                  Show this help message and exit
help
}

#parse parameteres
ARGS=""
re='^[0-9]+$'
THREADS="6"
OVERLAP="0.5"
RAM="8"
TE_length="80"
LENGTH="80"

while [[ $# -gt 0 ]]; do
	case $1 in
		--mate1)
		MATE1=$2
		shift 2
		;;
		--mate2)
		MATE2=$2
		shift 2
		;;
		--stranded)
		STRANDED=$2
					if [[ -z "$STRANDED" ]]; then
							echo "$(tput setaf 1)WARNING: The strand-specific parameter was not selected! Please use -std | --stranded (RF or FR)	Exiting..." && exit 1; else
							if [[ $STRANDED != "RF" ]] && [[ $STRANDED != "FR" ]]; then
								echo -e "\nERROR: -$STRANDED- is not a valid stranded option. Exiting..." && exit 1
							fi
					fi
		shift 2
		;;
		--transcripts)
		TRANSCRIPTS=$2
		shift 2
		;;
		--output)
		OUTPUT=$2
						if [[ -z "$OUTPUT" ]]; then
								OUTPUT=$(pwd)
						fi
						if [[ "$OUTPUT" == */ ]]; then
							OUTPUT=$(sed 's/.$//g' <<< "$OUTPUT")
						fi
		shift 2
		;;
		--ref_TEs)
		REF_TEs=$2
		shift 2
		;;
		--threads)
		THREADS=$2
            if ! [[ $THREADS =~ $re ]] ; then
                echo "ERROR: -$2- is not a valid value to threads. Exiting..."
                exit 1
            fi
		shift 2
		;;
		--ram)
		RAM=$2
						if ! [[ $RAM =~ $re ]] ; then
								echo "ERROR: $2 is not a valid value to ram. Exiting..."
								exit 1
						fi
		shift 2
		;;
		--TE_length)
		TE_length=$2
						if ! [[ $TE_length =~ $re ]] ; then
								echo "ERROR: $2 is not a valid value to --TE_length. Exiting..."
								exit 1
						fi
		shift 2
		;;
		--min_length)
		LENGTH=$2
						if ! [[ $LENGTH =~ $re ]] ; then
								echo "ERROR: $2 is not a valid value to --length. Exiting..."
								exit 1
						fi
		shift 2
		;;
		--overlap)
		OVERLAP=$2
		shift 2
		;;
		--min_TPM)
		TPM=$2
						if ! [[ $TPM =~ $re ]] ; then
								echo "ERROR: $2 is not a valid value to --min_TPM. Exiting..."
								exit 1
						fi
						if [[ -z "$TPM" ]]; then
								TPM="1"
						fi
		shift 2
		;;
		-h | --help)
			usage
			exit 1
			;;
		-*|--*)
			echo "Unknown option $1"
			exit 1
			;;
		*)
			ARGS+=("$1")
			shift
			;;
	esac
done

setting_parameters () {
echo "mate1 = $MATE1"
echo "mate2 = $MATE2"
echo "stranded = $STRANDED"
echo "out dir = $OUTPUT"
echo "TEs ref = $REF_TEs"
echo "threads = $THREADS"
echo "RAM = $RAM"
echo "OVERLAP = $OVERLAP"

}
setting_parameters

mkdir "$OUTPUT"/tmp 2>/dev/null
r1_samples=$(sed 's/,/\n/g' <<< "$MATE1")
r2_samples=$(sed 's/,/\n/g' <<< "$MATE2")
paste <(echo "$r1_samples") <(echo "$r2_samples") --delimiters '\t' > "$OUTPUT"/tmp/samples.lst

while IFS= read -r replicates; do
	mate1_rep=$(cut -f1 <<< "$replicates")
	mate2_rep=$(cut -f2 <<< "$replicates")
	sample=$(cut -f1 <<< "$replicates" | sed 's|.*/||g; s|/.*/||g; s/.fastq//g; s/.fq//g; s/.gz//g; s/_R1//g')
	echo "$sample"
	mkdir "$OUTPUT"/"$sample" 2>/dev/null
	#./scripts/part1-mode3.sh
	#./scripts/part2-mode3.sh
	./scripts/part3-mode3.sh
done < "$OUTPUT"/tmp/samples.lst

############################################

chimeric_list () {
echo "Chimeric reads coverage..."

isoforms=$(cut -f1 "$OUTPUT"/all_candidates.tsv | sort | uniq)

while IFS= read -r IDs; do
	reads=$(grep "$IDs" "$OUTPUT"/all_candidates.tsv | cut -f2)
	while IFS= read -r id; do
		mate=$(grep '/1' <<< "$id")
			if [ ! -z "$mate" ]; then
				oppsit=$(sed 's|/1|/2|g' <<< "$mate")
				chim_read=$(grep "$oppsit" "$OUTPUT"/transcripts_chim.lst | grep -w "$IDs")
				else
				mate=$(grep '/2' <<< "$id")
				oppsit=$(sed 's|/2|/1|g' <<< "$mate")
				chim_read=$(grep "$oppsit" "$OUTPUT"/transcripts_chim.lst | grep -w "$IDs")
			fi
			if [ ! -z "$chim_read" ]; then
				TE_chim=$(grep -w "$id" "$OUTPUT"/all_candidates.tsv | cut -f4)
				isoform=$(cut -f1 <<< "$chim_read")
				printf "$isoform""\t""$TE_chim""\n" >> "$OUTPUT"/coverage.rm
				chim_read=""
			fi
	done <<< "$reads"

	if [ -s "$OUTPUT"/coverage.rm ]; then
		chim_isoform=$(cut -f1 "$OUTPUT"/coverage.rm | head -1)
		TE=$(cut -f2 coverage.rm | sort | uniq -c | sort -V -k 1,1 | awk '{print $2,$1}' | sed 's/ /\t/g' | head -1)
		printf "$chim_isoform""\t""$TE""\n" >> "$OUTPUT"/chimeric_transcripts.ct
		rm "$OUTPUT"/coverage.rm
	fi

done <<< "$isoforms"

echo "done"
}

gene_identification () {

cut -f1 "$OUTPUT"/chimeric_transcripts.ct > "$OUTPUT"/chim_trinity_IDs.lst
seqtk "$OUTPUT"/subseq "$OUTPUT"/trinity_out_dir/Trinity.fasta "$OUTPUT"/chim_trinity_IDs.lst > "$OUTPUT"/chim_trinity.fa

makeblastdb -in "$TRANSCRIPTS" -dbtype nucl -out "$OUTPUT"/transcripts_db
blastn -query "$OUTPUT"/chim_trinity.fa -db "$OUTPUT"/transcripts_db -outfmt "6 qseqid sseqid length pident gaps mismatch qlen slen qstart qend sstart send evalue bitscore" -num_threads "$THREADS" > "$OUTPUT"/blast-result.tsv

isoform_IDs=$(cut -f1 "$OUTPUT"/blast-result.tsv | sort | uniq)

while IFS= read -r line; do
        isoform_length=$(grep "$line" "$OUTPUT"/blast-result.tsv | cut -f7 | head -1)
        best_hit_ID=$(grep "$line" "$OUTPUT"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f2 | head -1)
        best_hit_length=$(grep "$line" "$OUTPUT"/blast-result.tsv | grep "$best_hit_ID" | awk '{sum+=$3;} END{print sum;}')
	identity=$(grep "$line" "$OUTPUT"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f4 | head -1)
	ref_len=$(grep "$line" "$OUTPUT"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f8 | head -1)
	chim_len=$(grep "$line" "$OUTPUT"/blast-result.tsv | sort -k1,1 -k14,14gr | cut -f7 | head -1)
	transcript_gene=$(sed 's/_/\t/g' <<< "$best_hit_ID")
        perc=$(( best_hit_length*100/isoform_length ))
                if [ "$perc" -gt "$LENGTH" ]; then
                        printf "$line""\n" >> "$OUTPUT"/IDs_isoforms.lst
                        printf "$best_hit_ID""\n" >> "$OUTPUT"/IDs_genes.lst
			printf "$line""\t""$transcript_gene""\t""$identity""\t""$chim_len""\t""$ref_len""\t""$best_hit_length""\n" >> "$OUTPUT"/blast_matches.tsv
                fi
done <<< "$isoform_IDs"

paste "$OUTPUT"/IDs_isoforms.lst "$OUTPUT"/IDs_genes.lst > "$OUTPUT"/chim_IDs.lst
sed -i 's/^/>/g' "$OUTPUT"/IDs_genes.lst

seqtk subseq "$OUTPUT"/chim_trinity.fa "$OUTPUT"/IDs_isoforms.lst > "$OUTPUT"/gene_homology.fa
}

function renaming() {

python3 - << END
fasta= open('"$OUTPUT"/gene_homology.fa')
newnames= open('"$OUTPUT"/IDs_genes.lst')
newfasta= open('"$OUTPUT"/gene_homology-newIDs.fa', 'w')

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

rm IDs_isoforms.lst IDs_genes.lst
}

output () {
isoforms=$(cut -f1 chimeric_transcripts.ct)

while IFS= read -r line; do
	TE_cov=$(grep -w "$line" chimeric_transcripts.ct | cut -f2,3)
	blast_result=$(grep -w "$line" blast_matches.tsv)
	if [[ ! -z "$blast_result" ]]; then
		printf "$blast_result""\t""$TE_cov""\n" >> output.tsv
	fi
done <<< "$isoforms"

}
#massive_data
#managing_files
#chimeric_list
#gene_identification
#renaming
#output

#cut -f3,8 output.tsv | sed 's/\t/_/g; s/_I//g; s/-LTR_DM//g; s/_DM//g; s/-I//g; s/_LTR//g'
