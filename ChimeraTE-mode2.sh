#!/bin/bash

set -a

function usage() {

bash ./scripts/mode2/header.sh
   cat << help

USAGE:

ChimTE-mode2.sh [--mate1 <mate1_replicate1.fastq.gz,mate1_replicate2.fastq.gz>] [--mate2 <mate2_replicate1.fastq.gz,mate2_replicate2.fastq.gz>] [--te <TE_insertions.fa>] [--transcripts <transcripts.fa>] [--stranded <rf-stranded or fr-stranded>] [--project <project_name>] [options]

#Mandatory arguments:

  --mate1                 mate 1 from paired-end reads

  --mate2                 mate 2 from paired-end reads

  --te                    TE insertions (fasta)

  --transcripts           transcripts (fasta)

  --stranded              Select "rf-stranded" if your reads are reverse->forward; or "fr-stranded" if they are forward->reverse

  --project               project name

#Optional arguments:

  --fpkm                  Minimum FPKM expression (default = 1)

  --coverage              Minimum chimeric reads as support (default = 2)

  --threads               Number of threads, (default = 6)

  --assembly              Perform transcripts assembly -HIGH time consuming- (default = deactivated)
       |
       |
       -------> Mandatory if --assembly is activated:
       |        --ref_TEs             "species" database used by RepeatMasker (flies, human, mouse, arabidopsis);
       |                              or a built TE library in fasta format
       |
       |
       -------> Optional if --assembly is activated:
                --ram                 RAM memory (default: 8 Gbytes)
                --overlap             Minmum overlap (0.1 to 1) between read length and TE insertion (default = 0.5)
                --TE_length           Minimum TE length to keep it from RepeatMasker output (default = 80bp)
                --min_length          Minimum identity between de novo assembled transcripts and reference transcripts (default = 80%)

help
}

#parse parameteres
ARGS=""
RED='\033[1;31m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
NC='\033[0m'
re='^[0-9]+$'
COV="2"
THREADS="6"
OVERLAP="0.5"
RAM="8"
TE_length="80"
LENGTH="80"
DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
#folder=$(date | awk '{print $2,"-",$3,"-",$4,"-",$5}' | sed 's/ //g; s/:/\t/g' | awk '{print $1,$2}' | sed 's/ /h/g')

bash ./scripts/mode2/header.sh

if [ "$#" -eq 0 ]; then echo -e "\n${RED}ERROR${NC}: No parameters provided! Exiting..." && usage >&2; exit 1; fi

while [[ $# -gt 0 ]]; do
	case $1 in
		--mate1)
		MATE1=$2
		shift 2;;
		--mate2)
		MATE2=$2
		shift 2;;
    --te)
    TE_FASTA=$2
    shift 2;;
    --transcripts)
    TRANSCRIPTS_FASTA=$2
    shift 2;;
    --stranded)
    STRANDED=$2
    shift 2;;
    --project)
    PROJECT="$DIR"/projects/$2
    shift 2;;
    --assembly)
    ASSEMBLY=$#
    shift 1;;
    --ram)
		RAM=$2
    if ! [[ $RAM =~ $re ]] ; then
      echo -e "\n${RED}ERROR${NC}: --ram "$RAM" is not a valid value!! Use >= 8. Exiting..." &&  exit 1; else
      if [[ $RAM -lt "8" ]] ; then
        echo -e "\n${RED}ERROR${NC}: --ram "$RAM" is not a valid value!! Use >= 8. Exiting..." &&  exit 1
      fi
    fi
		shift 2;;
		--TE_length)
		TE_length=$2
			if ! [[ $TE_length =~ $re ]] ; then
				echo -e "\n${RED}ERROR${NC}: $TE_length is not a valid value to --TE_length. Exiting..." &&  exit 1
			fi
		shift 2;;
		--min_length)
		LENGTH=$2
    if ! [[ $LENGTH =~ $re ]] ; then
      echo -e "\n${RED}ERROR${NC}: $LENGTH is not a valid value to --TE_length. Exiting..." &&  exit 1
    fi
		shift 2;;
    --fpkm)
		FPKM=$2
    if ! [[ $FPKM =~ $re ]] ; then
      echo -e "\n${RED}ERROR${NC}: --fpkm "$FPKM" is not a valid value!! Use >= 2. Exiting..." &&  exit 1; else
      if [[ $FPKM -lt "1" ]] ; then
        echo -e "\n${RED}ERROR${NC}: --fpkm "$FPKM" is not a valid value!! Use >= 2. Exiting..." &&  exit 1
      fi
    fi
		shift 2;;
    --ref_TEs)
    REF_TEs=$2
    shift 2;;
    --coverage)
    COV=$2
    if ! [[ $COV =~ $re ]] ; then
      echo -e "\n${RED}ERROR${NC}: --coverage "$COV" is not a valid value!! Use >= 2. Exiting..." &&  exit 1; else
      if [[ $COV -lt "2" ]] ; then
        echo -e "\n${RED}ERROR${NC}: --coverage "$COV" is not a valid value!! Use >= 2. Exiting..." &&  exit 1
      fi
    fi
    shift 2;;
    --overlap)
    OVERLAP=$2
    if ! [[ $OVERLAP =~ ^[0-9]+(\.[0-9]+)?$ ]]; then #&& [ "$OVERLAP" =~ $re ]; then
      echo -e "\n${RED}ERROR${NC}: --coverage "$OVERLAP" is not a valid value!! Use a value between 0.1 to 1. Exiting..." &&  exit 1; else
      if [[ $OVERLAP -gt "1" ]] ; then
        echo -e "\n${RED}ERROR${NC}: --coverage "$OVERLAP" is not a valid value -if2!! Use a value between 0.1 to 1. Exiting..." &&  exit 1
      fi
    fi
    shift 2;;
    --threads)
    THREADS=$2
    if ! [[ $THREADS =~ $re ]] ; then
      echo -e "\n${RED}ERROR${NC}: --threads "$THREADS" is not a valid value!! Use >= 6. Exiting..." &&  exit 1; else
      if [[ $THREADS -lt "6" ]] ; then
        echo -e "\n${RED}ERROR${NC}: --threads "$THREADS" is not a valid value!! Use >= 6. Exiting..." && exit 1
      fi
    fi
    shift 2;;
    -h | --help)
			usage
			exit 1;;
		-*|--*)
			echo "Unknown option $1"
			exit 1;;
		*)
			ARGS+=("$1")
      usage
			shift;;

	esac
done

echo -ne "\nChecking files..."

if [[ -z "$MATE1" ]]; then
  echo -e "\n${RED}ERROR${NC}: --mate1 fastq files were not provided! Exiting..." && exit 0; else
    r1_samples=$(sed 's/,/\n/g' <<< "$MATE1")
    while IFS= read -r mate; do
      if [[ ! -f "$mate" ]]; then
        echo -e "\n${RED}ERROR${NC}: $mate provided with --mate1 was not found!! Exiting..." && exit 1; else
        if !(egrep -q '.fastq|.fq') <<< "$mate"; then
          echo -e "\n${RED}ERROR${NC}: Files provided with --mate1 are not in .fastq or .fq extension! Exiting..." && exit 1
        fi
      fi
    done <<< "$r1_samples"
fi

if [[ -z "$MATE2" ]]; then
  echo -e "\n${RED}ERROR${NC}: --mate2 fastq files were not provided! Exiting..." && exit 1; else
  r2_samples=$(sed 's/,/\n/g' <<< "$MATE2")
  while IFS= read -r mate; do
    if [[ ! -f "$mate" ]]; then
      echo -e "\n${RED}ERROR${NC}: $mate provided with --mate2 was not found!! Exiting..." && exit 1; else
      if !(egrep -q '.fastq|.fq') <<< "$mate"; then
        echo -e "\n${RED}ERROR${NC}: Files provided with --mate2 are not in .fastq or .fq extension! Exiting..." && exit 1
      fi
    fi
  done <<< "$r2_samples"
fi

if [[ -z "$TE_FASTA" ]]; then
  echo -e "\n${RED}ERROR${NC}: --te fasta file was not provided! Exiting..." && exit 1; else
  if [[ ! -f "$TE_FASTA" ]]; then
    echo -e "\n${RED}ERROR${NC}: $TE_FASTA provided with --te was not found!! Exiting..." && exit 1; else
    if !(egrep -q '.fasta|.fa') <<< "$TE_FASTA"; then
      echo -e "\n${RED}ERROR${NC}: $TE_FASTA provided with --te is not in .fasta extension! Exiting..." && exit 1
    fi
  fi
fi

if [[ -z "$TRANSCRIPTS_FASTA" ]]; then
  echo -e "\n${RED}ERROR${NC}: --transcripts fasta file was not provided! Exiting..." && exit 1; else
  if [[ ! -f "$TRANSCRIPTS_FASTA" ]]; then
    echo -e "\n${RED}ERROR${NC}: $TRANSCRIPTS_FASTA provided with --transcripts was not found!! Exiting..." && exit 1; else
    if !(egrep -q '.fasta|.fa') <<< "$TRANSCRIPTS_FASTA"; then
      echo -e "\n${RED}ERROR${NC}: $TRANSCRIPTS_FASTA provided with --transcripts is not in .fasta extension! Exiting..." && exit 1
    fi
  fi
fi

if [[ -z "$STRANDED" ]]; then
        echo -e "\n${RED}ERROR${NC}: The strand-specific parameter was not selected! Please use -s | --strand (rf-stranded or fr-stranded)	Exiting..." && exit 1; else
        if [ "$STRANDED" != "rf-stranded" ] && [ "$STRANDED" != "fr-stranded" ]; then
                echo -e "\n${RED}ERROR${NC}: The option -$STRANDED- is not accepted, please use --strand (rf-stranded or fr-stranded)	Exiting..." && exit 1
        fi
fi

echo -e "\t${GREEN}==> DONE!${NC}"

echo -e "\nParameters:
==> mate1 = $MATE1
==> mate2 = $MATE2
==> TEs fasta = $TE_FASTA
==> Transcripts fasta = $TRANSCRIPTS_FASTA
==> Project = $PROJECT
==> strand = $STRANDED
==> Coverage = $COV"

if [ ! -z $ASSEMBLY ]; then
  echo -e "==> assembly = ON\n==> RAM = $RAM\n==> TE_length = $TE_length\n==> min_length = $LENGTH"
fi

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

	#echo -e "\nRunning part1.sh - Gene and TE alignments $sample"
	#./scripts/mode2/part1-mode2.sh

	#echo -e "\nRunning part2.sh - Identifying all chimeric read pairs $sample_names"
	#./scripts/mode2/part2-mode2.sh

  #echo -e "\nRunning expression.sh - calculating expression level of chimeric transcripts $sample_names"
  #./scripts/mode2/expression-mode2.sh

  if [[ -f "$PROJECT"/"$sample"/"$sample"_output-final.ct ]]; then
    cp "$PROJECT"/"$sample"/"$sample"_output-final.ct "$PROJECT"/tmp
  fi

  #if [ ! -z $ASSEMBLY ]; then
  #  mkdir "$PROJECT"/"$sample"/trinity_out; TRINITY_OUT="$PROJECT/$sample/trinity_out"
  #  ./scripts/mode2/assembly.sh
  #  ./scripts/mode2/assembled-chimeras.sh; else
  #  echo "--assembly option is OFF"
  #fi
done < "$PROJECT"/tmp/samples.lst
fi

#echo -e "\nRunning replicability.sh - Checking for chimeras found in more than 1 RNA-seq replicate"
#./scripts/mode2/replicability-mode2.sh
