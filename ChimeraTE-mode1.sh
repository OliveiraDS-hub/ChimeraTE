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

   --project	project name

#Optional arguments:

   --window	Upstream and downstream window size (default = 3000)

   --overlap  Minimum overlap between chimeric reads and TE insertions

   --utr 	It must be used if your gene annotation (-a | --gene) has UTR regions (default = off)

   --fpkm   Minimum fpkm to consider a gene as expressed (default = 1)

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
    PROJECT=projects/$2
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
echo "genome = $GENOME"
echo "TE gtf = $TE_ANNOTATION"
echo "gene gtf = $GENE_ANNOTATION"
echo "project = $PROJECT"
echo "window = $WINDOW"
echo "overlap = $OVERLAP"
echo "fpkm = $FPKM"

mkdir "$PROJECT"/tmp 2>/dev/null
r1_samples=$(sed 's/,/\n/g' <<< "$MATE1")
r2_samples=$(sed 's/,/\n/g' <<< "$MATE2")
paste <(echo "$r1_samples") <(echo "$r2_samples") --delimiters '\t' > "$PROJECT"/tmp/samples.lst


while IFS= read -r replicates; do
	mate1_rep=$(cut -f1 <<< "$replicates")
	mate2_rep=$(cut -f2 <<< "$replicates")
  sample=$(cut -f1 <<< "$replicates" | sed 's|.*/||g; s|/.*/||g; s/.fastq//g; s/.fq//g; s/.gz//g; s/_R1//g')
  echo "$sample"
  mkdir "$PROJECT" "$PROJECT"/"$sample" "$PROJECT"/"$sample"/TE_info "$PROJECT"/"$sample"/gene_info "$PROJECT"/"$sample"/alignment "$PROJECT"/"$sample"/tmp "$PROJECT"/"$sample"/gene_info/tmp "$PROJECT"/"$sample"/plot 2>/dev/null
  export TMP="$PROJECT"/"$sample"/gene_info/tmp
  export ALN="$PROJECT"/"$sample"/alignment
  export BAM_FILE="$PROJECT"/"$sample"/alignment/accepted_hits.bam
  export TE_info="$PROJECT"/"$sample"/TE_info
  export GENE_info="$PROJECT"/"$sample"/gene_info

  #echo -e "\nRunning strd_specific.sh - Strand-specific alignment management"
  #./scripts/mode1/strd_specific.sh

  #echo -e "\nRunning exp_TEs.sh - Identifying expressed TE insertions -\n"
  #./scripts/mode1/exp_TEs.sh

  echo -e "\nRunning te-initiated-search.sh  - Searching for TE-initiated transcripts -\n"
  ./scripts/mode1/te-initiated-search.sh
  #echo -e "TE-initiated searching has finished sucessfully -- check TE-initiated.ct file\n"

  #echo -e "\nRunning te-terminated-search.sh  - Searching for TE-terminated transcripts -\n"
  #./scripts/mode1/te-terminated-search.sh
  #echo -e "\nTE-terminated searching has finished sucessfully -- check TE-terminated.ct file\n"

  #echo -e "\nRunning te-exonized-search.sh  - Searching for TE-exonized transcripts -\n"
  #./scripts/mode1/te-exonized-search.sh
  #echo -e "\nTE-exonized searching has finished sucessfully -- check TE-exonized.ct file\n"

  #echo -e "\n- Managing output tables -\n"
  #./scripts/mode1/te-init-out.sh
  #./scripts/mode1/te-exon-out.sh
  #./scripts/mode1/te-term-out.sh

  #echo -e "Plotting... \n"
  #./scripts/mode1/plot.sh

  #echo -e "\nThe process has finished successfully"
done < "$PROJECT"/tmp/samples.lst
