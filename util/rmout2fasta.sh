#!/bin/bash

function usage() {
   cat << help

Conversion of .out table from RepeatMasker to fasta file used by ChimeraTE - mode 2

#Mandatory arguments:

  --genome    file with genome (.fasta)
  --rm   file from RepeatMasker (.out)
  --out   output with TE insertions (.fasta)
help
}

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

egrep -v 'Satellite|Simple_repeat|rRNA|Low_complexity|RNA|ARTEFACT' "$RM" | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | tail -n +4 | awk '{Sense=$9;sub(/C/,"-",Sense);$9=Sense; Family=$10;gsub(/_D.*|-int|_I|-I|_LTR|-LTR/,"",Family);$10=Family ;print $5"\t"$6"\t"$7"\t"$10"\t"$1"\t"$9}' | sed 's/ /\t/g' > rmasker.bed
bedtools getfasta -fi "$GENOME" -bed rmasker.bed -s -name | sed 's/(.*)//g; s/::.*//g' > "$OUTPUT"
rm rmasker.bed *fai
