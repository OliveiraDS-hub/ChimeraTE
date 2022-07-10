#!/bin/bash

set -a

function usage() {
   cat << help

######################################
######## ChimeraTE v1.0 mode2 ########
######################################

USAGE:

-One-replicate:

ChimTE-mode2.sh [--mate1 <mate1.fastq.gz>] [--mate2 <mate2.fastq.gz>] [--te <TE_insertions.fa>] [--transcripts <transcripts.fa>] [--stranded <rf-stranded or fr-stranded>] [--project <project_name>] [options]

-Multi-replicates:

ChimTE-mode2.sh [--mate1 <mate1_replicate1.fastq.gz,mate1_replicate2.fastq.gz>] [--mate2 <mate2_replicate1.fastq.gz,mate2_replicate2.fastq.gz>] [--te <TE_insertions.fa>] [--transcripts <transcripts.fa>] [--stranded <rf-stranded or fr-stranded>] [--project <project_name>] [options]


#Mandatory arguments:

  --mate1 mate 1 from paired-end reads

  --mate2 mate 2 from paired-end reads

  --te  TE insertions (fasta)

  --transcripts transcripts (fasta)

  --stranded  Select "rf-stranded" if your reads are reverse->forward; or "fr-stranded" if they are forward->reverse

  --project project name

#Optional arguments:

  --coverage Minimum chimeric reads as support (default: 2)

  --threads Number of threads, (default: 6)

help
}

#parse parameteres
ARGS=""
re='^[0-9]+$'
COV="2"
THREADS="6"
DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

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
    --te)
    TE_FASTA=$2
    shift 2
    ;;
    --transcripts)
    TRANSCRIPTS_FASTA=$2
    shift 2
    ;;
    --stranded)
    STRANDED=$2
    if [ -z "$STRANDED" ]; then
            echo "WARNING: The strand-specific direction was not selected! Please use -s | --strand (rf-stranded or fr-stranded)	Exiting..." && exit 1; else
            if [ "$STRANDED" == "rf-stranded" ]; then
                    echo "The $STRANDED accept only reverse->forward alignments"; else
                    if [ "$STRANDED" == "fr-stranded" ]; then
                            echo "The $STRANDED accept only reverse->forward alignments"; else
                            echo "WARNING: $STRANDED is not an option to --strand parameter! Exiting..." && exit 1
                    fi
            fi
    fi
    shift 2
    ;;
    --project)
    PROJECT="$DIR"/projects/$2
    shift 2
    ;;
    --coverage)
    COV=$2
    shift 2
    ;;
    --threads)
    THREADS=$2
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

echo "mate1 = $MATE1"
echo "mate2 = $MATE2"
echo "TEs fasta = $TE_FASTA"
echo "Transcripts fasta = $TRANSCRIPTS_FASTA"
echo "Project = $PROJECT"
echo "strand = $STRANDED"
echo "Coverage = $COV"

if [[ ! -d "$PROJECT" || ! -d "$PROJECT"/tmp ]]; then
  mkdir "$PROJECT" 2>/dev/null
  mkdir "$PROJECT"/tmp 2>/dev/null
fi

r1_samples=$(sed 's/,/\n/g' <<< "$MATE1")
r2_samples=$(sed 's/,/\n/g' <<< "$MATE2")
paste <(echo "$r1_samples") <(echo "$r2_samples") --delimiters '\t' > "$PROJECT"/tmp/samples.lst

if [[ ! -z $r1_samples ]]; then

while IFS= read -r SAMP; do
	mate1=$(cut -f1 <<< "$SAMP")
	mate2=$(cut -f2 <<< "$SAMP")
	sample=$(cut -f1 <<< "$SAMP" | sed 's|.*/||g; s|/.*/||g; s/.fastq//g; s/.fq//g; s/gz//g; s/_R1//g')

  if [[ ! -d "$PROJECT"/"$sample" ]]; then
	mkdir "$PROJECT"/"$sample" "$PROJECT"/indexes/ "$PROJECT"/"$sample"/alignment "$PROJECT"/"$sample"/chimeric_reads/ "$PROJECT"/"$sample"/plot "$PROJECT"/"$sample"/tmp_dir 2>/dev/null
  fi

	export ALN="$PROJECT"/"$sample"/alignment
	export READS="$PROJECT"/"$sample"/chimeric_reads

	echo -e "\nRunning part1.sh - Gene and TE alignments $sample"
	./scripts/mode2/part1-mode2.sh

	echo -e "\nRunning part2.sh - Identifying all chimeric read pairs $sample_names"
	./scripts/mode2/part2-mode2.sh

  echo -e "\nRunning expression.sh - calculating expression level of chimeric transcripts $sample_names"
  ./scripts/mode2/expression-mode2.sh

  if [[ -f "$PROJECT"/"$sample"/"$sample"_output-final.ct ]]; then
    cp "$PROJECT"/"$sample"/"$sample"_output-final.ct "$PROJECT"/tmp
  fi
done < "$PROJECT"/tmp/samples.lst
fi

echo -e "\nRunning replicability.sh - Checking for chimeras found in more than 1 RNA-seq replicate"
./scripts/mode2/replicability-mode2.sh
