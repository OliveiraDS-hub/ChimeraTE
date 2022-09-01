#!/bin/bash

function usage() {
   cat << help

Conversion of .out table from RepeatMasker to fasta file used by ChimeraTE - mode 2

#Mandatory arguments:

  --genome    	file with genome (.fasta)
  --rm   	file from RepeatMasker (.out)
  --out   	output with TE insertions (.fasta)
help
}

RED='\033[1;31m'
NC='\033[0m'
while [[ $# -gt 0 ]]; do
    case $1 in
    --genome)
    GENOME=$2
    shift 2
    ;;
    --rm)
    RM=$2
    shift 2
    ;;
    --out)
    OUTPUT=$2
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

if [[ -z "$GENOME" ]]; then
  echo -e "\n${RED}ERROR${NC}: --genome fasta file was not provided! Exiting..." && exit 1; else
  if [[ ! -f "$GENOME" ]]; then
    echo -e "\n${RED}ERROR${NC}: $GENOME provided with --genome does not exist!! Exiting..." && exit 1; else
    if !(egrep -q '.fasta|.fa|.fna') <<< "$GENOME"; then
      echo -e "\n${RED}ERROR${NC}: $GENOME provided with --genome is not in .fasta extension! Exiting..." && exit 1
    fi
  fi
fi

if [[ -z "$RM" ]]; then
  echo -e "\n${RED}ERROR${NC}: --rm .out file was not provided! Exiting..." && exit 1; else
  if [[ ! -f "$RM" ]]; then
    echo -e "\n${RED}ERROR${NC}: $RM provided with --rm does not exist!! Exiting..." && exit 1
  fi
fi

if [[ -z "$OUTPUT" ]]; then
  echo -e "\n${RED}ERROR${NC}: --out file name was not provided! Exiting..." && exit 1
fi


egrep -v 'Satellite|Simple_repeat|rRNA|Low_complexity|RNA|ARTEFACT' "$RM" | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | tail -n +4 | awk '{Sense=$9;sub(/C/,"-",Sense);$9=Sense; Family=$10;gsub(/_D.*|-int|_I|-I|_LTR|-LTR/,"",Family);$10=Family ;print $5"\t"$6"\t"$7"\t"$10"\t"$1"\t"$9}' | sed 's/ /\t/g' > rmasker.bed
bedtools getfasta -fi "$GENOME" -bed rmasker.bed -s -name | sed 's/(.*)//g; s/::.*//g' > "$OUTPUT"
rm rmasker.bed *fai
