#!/bin/bash
set -a

function usage() {
   cat << help
######################################
#####   ChimeraTE v1.0 - MODE 1  #####
######################################

USAGE:

-One-replicate:

ChimTE-mode1.sh [--mate1 <mate1.fastq.gz>] [--mate2 <mate2.fastq.gz>] [--genome <genome.fasta>] [--te <TE_insertions.gtf or TE_insertions.bed>] [--gene <gene_annotation.gtf>] [--project <project_name>]

#Mandatory arguments:

   --mate1	paired-end R1

   --mate2	paired-end R2

   --genome	genome sequence .fa

   --te	GTF file with TE coordinates

   --gene	GTF file with genes coordinates

   --project	project name (it's the name of the folder that will be created at the folder projects/)

#Optional arguments:

   --window	Upstream and downstream window size (default = 3000)

   --overlap  Minimum overlap between chimeric reads and TE insertions

   --utr 	It must be used if your gene annotation (-a | --gene) has UTR regions (default = off)

   --fpkm   Minimum fpkm to consider a gene as expressed (default = 1)

   --coverage	Minimum coverage for chimeric transcripts detection (default = 2)

   --replicate	Minimum recurrency of chimeric transcripts between RNA-seq replicates (default = 2)

   --threads  Number of threads (default:6)
help
}

#parse parameteres
ARGS=""
re='^[0-9]+$'
THREADS="6"
WINDOW="3000"
OVERLAP="0.5"
FPKM="1"
COV="2"
REP="2"
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
    --genome)
    GENOME=$2
    shift 2
    ;;
    --te)
    TE_ANNOTATION=$2
    shift 2
    ;;
    --gene)
    GENE_ANNOTATION=$2
    shift 2
    ;;
    --project)
    PROJECT="$DIR"/projects/$2
    shift 2
    ;;
    --window)
    WINDOW=$2
    shift 2
    ;;
    --overlap)
    OVERLAP=$2
    shift 2
    ;;
    --utr)
    UTR=$#
    shift 1
    ;;
    --fpkm)
    FPKM=$2
    shift 2
    ;;
    --coverage)
    COV=$2
    shift 2
    ;;
    --replicate)
    REP=$2
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

echo -e "\n######################################
#####   ChimeraTE v1.0 - MODE 1  #####
######################################

Parameters:

==> mate1:		$MATE1
==> mate2: 		$MATE2
==> genome: 		$GENOME
==> TE gtf: 		$TE_ANNOTATION
==> gene gtf: 		$GENE_ANNOTATION
==> project: 		$PROJECT
==> window: 		$WINDOW
==> overlap: 		$OVERLAP
==> fpkm:		$FPKM
==> replicability: 	$REP"

if [[ ! -d "$PROJECT" || ! -d "$PROJECT"/tmp ]]; then
  mkdir "$PROJECT" 2>/dev/null
  mkdir "$PROJECT"/tmp 2>/dev/null; else
  echo "The project $PROJECT has been found and it won't be overwritten"
fi

r1_samples=$(sed 's/,/\n/g' <<< "$MATE1")
r2_samples=$(sed 's/,/\n/g' <<< "$MATE2")
paste <(echo "$r1_samples") <(echo "$r2_samples") --delimiters '\t' > "$PROJECT"/tmp/samples.lst

while IFS= read -r replicates; do
	mate1_rep=$(cut -f1 <<< "$replicates")
	mate2_rep=$(cut -f2 <<< "$replicates")
  sample=$(cut -f1 <<< "$replicates" | sed 's|.*/||g; s|/.*/||g; s/.fastq//g; s/.fq//g; s/.gz//g; s/_R1//g')
  mkdir "$PROJECT" "$PROJECT"/"$sample" "$PROJECT"/"$sample"/TE_info "$PROJECT"/"$sample"/gene_info "$PROJECT"/"$sample"/alignment "$PROJECT"/"$sample"/tmp "$PROJECT"/"$sample"/gene_info/tmp "$PROJECT"/"$sample"/plot 2>/dev/null
  export TMP="$PROJECT"/"$sample"/gene_info/tmp
  export ALN="$PROJECT"/"$sample"/alignment
  export BAM_FILE="$PROJECT"/"$sample"/alignment/accepted_hits.bam
  export TE_info="$PROJECT"/"$sample"/TE_info
  export GENE_info="$PROJECT"/"$sample"/gene_info

  echo -e "\nRunning strd_specific.sh - Strand-specific alignment management"
  ./scripts/mode1/strd_specific.sh

  echo -e "\nRunning exp_TEs.sh - Identifying expressed TE insertions -\n"
  ./scripts/mode1/exp_TEs.sh

  echo -e "\nRunning te-initiated-search.sh  - Searching for TE-initiated transcripts -\n"
  ./scripts/mode1/te-initiated-search.sh
  
  echo -e "\nRunning te-terminated-search.sh  - Searching for TE-terminated transcripts -\n"
  ./scripts/mode1/te-terminated-search.sh

  echo -e "\nRunning te-exonized-search.sh  - Searching for TE-exonized transcripts -\n"
  ./scripts/mode1/te-exonized-search.sh

  #echo -e "Plotting... \n"
  #./scripts/mode1/plot.sh

done < "$PROJECT"/tmp/samples.lst

echo -e "\nRunning replicability-mode1.sh  - Searching for chimeras identified in at least $REP replicates"
./scripts/mode1/replicability-mode1.sh

echo -e "\nThe process has finished successfully"
