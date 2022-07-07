#!/bin/bash

set -e

function usage() {
   cat << help

Repeatmasker analysis and parsing with One Code to Find Them All:

#Mandatory arguments:

  --genome    file with genome (.fasta)

  --ref_TEs   "species" database used by RepeatMasker (flies, human, mouse, arabidopsis...)

  --out       output file

#Optional arguments:

  --threads   Number of threads, (default: 6)

  --dist      Distance in nt between LTR and Internal regions, as well as fragments from the
              same TE family that will be merged, (default: 50)

help
}

THREADS="6"
DIST="50"

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

declare -a requirements=("RepeatMasker" "octa/build_dictionary.pl" "octa/one_code_to_find_them_all.pl")
for req in "${requirements[@]}"; do
	command -v "$req" >> /dev/null 2>&1
	if [[ $? -ne 0 ]] ; then
	    echo -e "$req not found. Exiting...\n"
	    exit 1
	fi
done

echo -e "Running RepeatMasker... it may take long\n"

RepeatMasker "$GENOME" -species "$SPECIES" -cutoff 225 -nolow -norna -a -s -par "$THREADS" >> /dev/null 2>&1

if [[ -f "$GENOME".out ]]; then
  echo -e "RepeatMasker process has been finished!! output: "$GENOME".out\n"
  sleep 2
fi

echo -e "Running One Code to Find them all...\n"

perl octa/build_dictionary.pl --rm "$GENOME".out > "$OUTPUT".dictio >> /dev/null 2>&1
perl octa/one_code_to_find_them_all.pl --rm "$GENOME".out --ltr "$OUTPUT".dictio --unknown --fasta "$GENOME" --insert "$DIST" >> /dev/null 2>&1

mkdir tmp_data tmp_data/scaf_ltr tmp_data/scaf_transposons tmp_data/scaf_copynumber tmp_data/scaf_TEs >> /dev/null 2>&1

mv "$GENOME"*sons.csv tmp_data/scaf_transposons/ ; mv "$GENOME"*ltr.csv tmp_data/scaf_ltr/ ; mv "$GENOME"*number.csv tmp_data/scaf_copynumber/ ; mv "$GENOME"*sorted.csv tmp_data/scaf_TEs/

sed -i '/^$/d' tmp_data/scaf_TEs/*sorted.csv
sed -i '1d' tmp_data/scaf_TEs/*sorted.csv
cat tmp_data/scaf_TEs/*csv | grep '###' > "$OUTPUT".octa

awk '{Sense=$9;sub(/C/,"-",Sense);$9=Sense; print $5"\tRepeatMasker\tsimilarity\t"$6"\t"$7"\t"$2"\t"$9"\t.\t"$10}' "$OUTPUT".octa > "$OUTPUT"_RM_final.gtf
rm -R tmp_data "$GENOME"*length "$GENOME"*masked "$GENOME"*.align "$GENOME"*cat.gz *dictio "$OUTPUT".octa "$GENOME"*.log.txt >> /dev/null 2>&1

if [[ -f "$OUTPUT"_RM_final.gtf ]]; then
  echo -e "One Code to Find Them All has been finished!! output: "$OUTPUT"_RM_final.gtf\n"
fi

