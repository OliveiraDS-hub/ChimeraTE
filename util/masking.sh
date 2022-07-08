#!/bin/bash

set -e

function usage() {
   cat << help

This script is used to perform Repeatmasker analysis and parsing with One Code to Find Them All, 
providing the gtf file with TE insertions in a proper format to run ChimeraTE mode1 (genome-guided).

#Mandatory arguments:

  --genome    file with genome (.fasta)

  --ref_TEs   species database used by RepeatMasker ("flies", "human", "mouse", "arabidopsis"...); 
  	      or a file with your custom TE library (.fasta)

#Optional arguments:

  --out       output directory (if not provided, the output will be created in the current folder)

  --threads   Number of threads, (default = 6)

  --dist      Distance in nt between LTR and Internal regions, as well as fragments from the
              same TE family that will be merged by One Code to Find them all, (default = 50)

#Citations:

Bailly-Bechet, Marc, Annabelle Haudry, and Emmanuelle Lerat. "“One code to find them all”: a perl tool to conveniently parse RepeatMasker output files." Mobile DNA 5.1 (2014): 1-15.

Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015 <http://www.repeatmasker.org>. 

help
}

THREADS="6"
DIST="50"
OUTPUT=$(pwd)
DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

while [[ $# -gt 0 ]]; do
	case $1 in
	--genome)
	GENOME=$2
	shift 2
	;;
    	--ref_TEs)
    	SPECIES=$2
    	shift 2
    	;;
    	--out)
    	OUTPUT=$2
    	shift 2
    	;;
    	--threads)
    	THREADS=$2
    	shift 2
    	;;
    	--dist)
    	DIST=$2
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

echo -e "\nChimeraTE - masking.sh script\n"

FILE=$(basename "$GENOME")
OUTfile=$(basename "$GENOME" | sed 's/.fasta//g; s/.fa//g')

declare -a requirements=("RepeatMasker" ""$DIR"/octa/build_dictionary.pl" ""$DIR"/octa/one_code_to_find_them_all.pl")
for req in "${requirements[@]}"; do
	command -v "$req" >> /dev/null 2>&1
	if [[ $? -ne 0 ]] ; then
	    echo -e "$req not found. Exiting...\n"
	    exit 1
	fi
done

echo -e "Running RepeatMasker... it may take long\n"

if egrep -q '.fasta|.fa' <<< "$SPECIES"; then
        RepeatMasker "$GENOME" -lib "$SPECIES" -cutoff 225 -nolow -norna -a -s -par "$THREADS" -dir "$OUTPUT" >> /dev/null 2>&1; else
        RepeatMasker "$GENOME" -species "$SPECIES" -cutoff 225 -nolow -norna -a -s -par "$THREADS" -dir "$OUTPUT" >> /dev/null 2>&1
fi

if [[ -f "$OUTPUT"/"$FILE".out ]]; then
	echo -e "RepeatMasker process has been finished!! output: "$OUTPUT"/"$FILE".out\n"
	sleep 2

	echo -e "Running One Code to Find them all...\n"
	perl "$DIR"/octa/build_dictionary.pl --rm "$OUTPUT"/"$FILE".out > "$OUTPUT"/"$FILE".dictio >> /dev/null 2>&1
	perl "$DIR"/octa/one_code_to_find_them_all.pl --rm "$OUTPUT"/"$FILE".out --ltr "$OUTPUT"/"$FILE".dictio --unknown --fasta "$GENOME" --insert "$DIST" >> /dev/null 2>&1
	mkdir "$OUTPUT"/tmp_data "$OUTPUT"/tmp_data/scaf_ltr "$OUTPUT"/tmp_data/scaf_transposons "$OUTPUT"/tmp_data/scaf_copynumber "$OUTPUT"/tmp_data/scaf_TEs >> /dev/null 2>&1
	mv "$OUTPUT"/*sons.csv "$OUTPUT"/tmp_data/scaf_transposons/ ; mv "$OUTPUT"/*ltr.csv "$OUTPUT"/tmp_data/scaf_ltr/ ; mv "$OUTPUT"/*number.csv "$OUTPUT"/tmp_data/scaf_copynumber/ ; mv "$OUTPUT"/*sorted.csv "$OUTPUT"/tmp_data/scaf_TEs/
	sed -i '/^$/d' "$OUTPUT"/tmp_data/scaf_TEs/*sorted.csv; sed -i '1d' "$OUTPUT"/tmp_data/scaf_TEs/*sorted.csv
	cat "$OUTPUT"/tmp_data/scaf_TEs/*csv | grep '###' > "$OUTPUT"/"$FILE".octa

	awk '{Sense=$9;sub(/C/,"-",Sense);$9=Sense; print $5"\tRepeatMasker\tsimilarity\t"$6"\t"$7"\t"$2"\t"$9"\t.\t"$10}' "$OUTPUT"/"$FILE".octa > "$OUTPUT"/"$OUTfile"_RM_final.gtf
	rm -R "$OUTPUT"/tmp_data "$OUTPUT"/"$FILE"*length "$OUTPUT"/*masked "$OUTPUT"/*.align "$OUTPUT"/*cat.gz "$OUTPUT"/*dictio "$OUTPUT"/"$FILE".octa "$OUTPUT"/*.log.txt "$OUTPUT"/*fasta >> /dev/null 2>&1
	if [[ -f "$OUTPUT"/"$OUTfile"_RM_final.gtf ]]; then
		echo -e "One Code to Find Them All has been finished!! output: "$OUTPUT"/"$OUTfile"_RM_final.gtf\n"
	fi
else
	echo "RepeatMasker output was not found!"
fi
