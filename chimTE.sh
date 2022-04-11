#!/bin/bash

set -e

OPTS=`getopt -o 1:2:t:g:p:s:ch --long mate1:,mate2:,te:,gene:,project:,strandness:,cutoff,help --name "$0" -- "$@"`

if [ $? != 0 ] ; then echo "Failed parsing options." >&2 ; exit 1 ; fi

eval set -- "$OPTS"

usage=("ChimeraTE - Usage commandline

        Required arguments:

        -1 | --mate1    	paired-end R1

        -2 | --mate2    	paired-end R2

        -g | --gene     	gene sequences .fa

        -t | --te       	TE sequences .fa

        -p | --project  	project name

	-s | --strandness	Select rf-stranded if your reads are reverse->forward; or fr-stranded if they are forward->reverse

        Optional arguments:

	-c | --cutoff   	Minimum chimeric pairs")

while true; do
  case "$1" in
    -1 | --mate1 ) r1="$2"; shift 2;;
    -2 | --mate2 ) r2="$2"; shift 2;;
    -t | --te ) input_TEs="$2"; shift 2;;
    -g | --gene ) input_genes="$2"; shift 2;;
    -p | --project ) proj="$2"; shift 2;;
    -s | --strandness ) strandness="$2"; shift 2;;
    -c | --cutoff ) cutoff="$2"; shift ;;
    -h | --help ) printf "$usage""\n"; exit 1;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

echo -n "Checking the required files and arguments..."

if [ -z "$strandness" ]; then
        echo "$(tput setaf 1)WARNING: The strand-specific direction was not selected! Please use -s | --strand (rf-stranded or fr-stranded)	Exiting..." && exit 1; else
        if [ "$strandness" == "rf-stranded" ]; then
                echo "The $strandness accept only reverse->forward alignments"; else
                if [ "$strandness" == "fr-stranded" ]; then
                        strand_info="The $strandness accept only reverse->forward alignments"; else
                        echo "$(tput setaf 1)WARNING: $strandness is not an option to --strand parameter! Exiting..." && exit 1
                fi
        fi
fi

if [ -z "$cutoff" ]; then
	cutoff="2" && cutoff_info="The minimum amount of chimeric fragments to support a chimeric transcript is 2 (default)" ; else
        cutoff_info="The minimum amount of chimeric fragments to support a chimeric transcript is $cutoff"
fi

echo "
##################################################################
##########						##########
##########	      ChimeraTE analysis v1.0		##########
##########                                              ##########
##################################################################
"

r1_samples=$(echo "$r1" | sed 's/,/\n/g')
r2_samples=$(echo "$r2" | sed 's/,/\n/g')

paste <(echo "$r1_samples") <(echo "$r2_samples") --delimiters '\t' > samples.lst

r1_samples_n=$(echo "$r1" | sed 's/,/\n/g' | wc -l)

echo "-Your RNA-seq data is composed by: $r1_samples_n replicates"
echo "-$cutoff_info"
echo "-$strand_info"

if [[ ! -z $r1_samples ]]; then
	
	while IFS= read -r SAMP; do
	mate1=$(cut -f1 <<< "$SAMP")
	mate2=$(cut -f2 <<< "$SAMP")

	sample=$(cut -f1 <<< "$SAMP" | sed 's|.*/||g; s|/.*/||g; s/.fastq//g; s/.fq//g; s/.gz//g; s/_R1//g')
	echo "-Starting analysis for: $sample sample"
	
	mkdir projects/"$sample" projects/"$sample"/alignment projects/"$sample"/ids_crossing/ projects/"$sample"/chimeric_reads/ projects/"$sample"/chimTE_out/ projects/"$sample"/chimTE_out/sg_sg projects/"$sample"/chimTE_out/sg_conc/ projects/"$sample"/chimTE_out/conc_conc/ projects/"$sample"/plot projects/"$sample"/tmp_dir projects/"$sample"/alignment/fpkm_counts/
	
		export ALN=projects/"$sample"/alignment
		export READS=projects/"$sample"/chimeric_reads
		export cutoff proj r1 r2 input_TEs input_genes r1_r2 mate1 mate2 strandness sample

	##### ALIGNMENT AND BAM CONVERSIONS #####
		date
		printf "\nRunning part1.sh - Gene and TE alignments $sample-\n"
		./scripts/part1.sh 

	##### Reads identification #####
		date
		printf "\nRunning part2.sh - Identifying all chimeric read pairs $sample_names-\n"
		./scripts/part2.sh 
	done < samples.lst
fi

echo "The chimeric transcripts analysis has been finished!! Check the results in the *final.ct files"







