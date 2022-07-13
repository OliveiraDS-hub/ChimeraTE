#!/bin/bash

set -e
RED='\033[1;31m'
NC='\033[0m'
dot='`'
function usage() {

  echo -e "\n${NC}   ________    _                          ${RED}____________${NC}"
  echo -e "${NC}  / ____/ /_  (_)___ ___  ___  _________ ${RED}/_  __/ ____/${NC}"
  echo -e "${NC} / /   / __ \/ / __ $dot __\/ _ \/ ___/ __ $dot${RED}// / / __/${NC}"
  echo -e "${NC}/ /___/ / / / / / / / / /  __/ /  / /_/ /${RED}/ / / /___${NC}"
  echo -e "${NC}\____/_/ /_/_/_/ /_/ /_/\___/_/   \__,_/${RED}/_/ /_____/${NC}"
  echo -e "${NC}-. .-.   .-. .-.   .-. .-.   .-. .-.   .${RED}-. .-.   .-. .${NC}"
  echo -e "${NC}||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|${RED}||\|||\ /|||\|${NC}"
  echo -e "${NC}|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\||${RED}|/ \|||\|||/ \ ${NC}"
  echo -e "${NC}~   $dot-~ $dot-$dot   $dot-~ $dot-$dot   $dot-~ $dot-~   $dot-~ $dot-${RED}$dot   $dot-~ $dot-$dot   ${NC}"

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

if [ "$#" -eq 0 ]; then echo -e "\n${RED}ERROR${NC}: No parameters provided! Exiting..." && usage >&2; exit 1; fi

RED='\033[1;31m'
GREEN='\033[1;32m'
NC='\033[0m'
THREADS="6"
DIST="50"
OUTPUT=$(pwd)
DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

while [[ $# -gt 0 ]]; do
	case $1 in
	--genome)
	GENOME=$2
	shift 2;;
  --ref_TEs)
  REF_TEs=$2
  shift 2;;
  --out)
  OUTPUT=$2
  shift 2;;
  --threads)
  THREADS=$2
    if [[ $THREADS -lt "6" ]] ; then
      echo -e "\n${RED}ERROR${NC}: --threads "$THREADS" is not a valid value!! Use >= 6. Exiting..." &&  exit 1
    fi
  shift 2;;
  --dist)
  DIST=$2
  shift 2;;
  -h | --help)
  usage
  exit 1;;
  -*|--*)
  echo "Unknown option $1"
  exit 1;;
  *)
  ARGS+=("$1")
  shift;;
  	esac
done

echo -ne "\nChecking input..."

if [[ -z "$GENOME" ]]; then
  echo -e "\n${RED}ERROR${NC}: --genome fasta file was not provided! Exiting..." && exit 1; else
  if [[ ! -f "$GENOME" ]]; then
    echo -e "\n${RED}ERROR${NC}: $GENOME provided with --genome was not found!! Exiting..." && exit 1; else
    if !(egrep -q '.fasta|.fa') <<< "$GENOME"; then
      echo -e "\n${RED}ERROR${NC}: $GENOME provided with --genome is not in .fasta extension! Exiting..." && exit 1
    fi
  fi
fi

if [[ -z "$REF_TEs" ]]; then
  echo -e "\n${RED}ERROR${NC}: Neither fasta file or RM_database were provided to --ref_TEs! Exiting..." && exit 1
fi

echo -e "\t${GREEN}==> DONE!${NC}"


FILE=$(basename "$GENOME")
OUTfile=$(basename "$GENOME" | sed 's/.fasta//g; s/.fa//g')

echo -ne "\nRunning RepeatMasker... it may take long"
if egrep -q '.fasta|.fa' <<< "$REF_TEs"; then
        RepeatMasker "$GENOME" -lib "$REF_TEs" -cutoff 225 -nolow -norna -a -s -par "$THREADS" -dir "$OUTPUT" > /dev/null; else
        RepeatMasker "$GENOME" -REF_TEs "$REF_TEs" -cutoff 225 -nolow -norna -a -s -par "$THREADS" -dir "$OUTPUT" > /dev/null
fi
echo -e "\t${GREEN}==> DONE!${NC}"

if [[ -f "$OUTPUT"/"$FILE".out ]]; then
	echo -e "RepeatMasker process has been finished!! output: "$OUTPUT"/"$FILE".out\n"
	sleep 2

	echo -ne "Running One Code to Find them all..."
	perl "$DIR"/octa/build_dictionary.pl --rm "$OUTPUT"/"$FILE".out > "$OUTPUT"/"$FILE".dictio >> /dev/null 2>&1
	perl "$DIR"/octa/one_code_to_find_them_all.pl --rm "$OUTPUT"/"$FILE".out --ltr "$OUTPUT"/"$FILE".dictio --unknown --fasta "$GENOME" --insert "$DIST" >> /dev/null 2>&1
	mkdir "$OUTPUT"/tmp_data "$OUTPUT"/tmp_data/scaf_ltr "$OUTPUT"/tmp_data/scaf_transposons "$OUTPUT"/tmp_data/scaf_copynumber "$OUTPUT"/tmp_data/scaf_TEs >> /dev/null 2>&1
	mv "$OUTPUT"/*sons.csv "$OUTPUT"/tmp_data/scaf_transposons/ ; mv "$OUTPUT"/*ltr.csv "$OUTPUT"/tmp_data/scaf_ltr/ ; mv "$OUTPUT"/*number.csv "$OUTPUT"/tmp_data/scaf_copynumber/ ; mv "$OUTPUT"/*sorted.csv "$OUTPUT"/tmp_data/scaf_TEs/
	sed -i '/^$/d' "$OUTPUT"/tmp_data/scaf_TEs/*sorted.csv; sed -i '1d' "$OUTPUT"/tmp_data/scaf_TEs/*sorted.csv
	cat "$OUTPUT"/tmp_data/scaf_TEs/*csv | grep '###' > "$OUTPUT"/"$FILE".octa

	awk '{Sense=$9;sub(/C/,"-",Sense);$9=Sense; print $5"\tRepeatMasker\tsimilarity\t"$6"\t"$7"\t"$2"\t"$9"\t.\t"$10}' "$OUTPUT"/"$FILE".octa > "$OUTPUT"/"$OUTfile"_RM_final.gtf
	rm -R "$OUTPUT"/tmp_data "$OUTPUT"/"$FILE"*length "$OUTPUT"/*masked "$OUTPUT"/*.align "$OUTPUT"/*cat.gz "$OUTPUT"/*dictio "$OUTPUT"/"$FILE".octa "$OUTPUT"/*.log.txt "$OUTPUT"/*fasta >> /dev/null 2>&1
	if [[ -f "$OUTPUT"/"$OUTfile"_RM_final.gtf ]]; then
    echo -e "\t${GREEN}==> DONE!${NC}"
		echo -e "One Code to Find Them All has been finished!! output: "$OUTPUT"/"$OUTfile"_RM_final.gtf\n"
	fi
else
	echo "RepeatMasker output was not found!"
fi
