#!/bin/bash
set -e

function usage() {
cat << help

Conversion of ENSEMBL native transcript IDs to the ChimeraTE Mode 2 format

#Mandatory arguments:

  --transcripts     transcripts downloaded from ENSEMBL (.fasta)
                    i.e.: Drosophila_melanogaster.BDGP6.32.cdna.all (D. melanogaster)
                          Arabidopsis_thaliana.TAIR10.cdna.all.fa (A. thaliana)

  --out   output file name (.fasta)
help
}

RED='\033[1;31m'
NC='\033[0m'

if [ "$#" -eq 0 ]; then echo -e "\n${RED}ERROR${NC}: No parameters provided! Exiting..." && usage >&2; exit 1; fi

while [[ $# -gt 0 ]]; do
  case $1 in
    --transcripts)
    TRANSCRIPTS=$2
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

if [[ -z "$TRANSCRIPTS" ]]; then
  echo -e "\n${RED}ERROR${NC}: --trancripts fasta file was not provided! Exiting..." && exit 1; else
  if [[ ! -f "$TRANSCRIPTS" ]]; then
    echo -e "\n${RED}ERROR${NC}: $TRANSCRIPTS provided with --transcripts does not exist!! Exiting..." && exit 1; else
    if !(egrep -q '.fasta|.fa|.fna') <<< "$TRANSCRIPTS"; then
      echo -e "\n${RED}ERROR${NC}: $TRANSCRIPTS provided with --transcripts is not in .fasta extension! Exiting..." && exit 1
    fi
  fi
fi

if [[ -z "$OUTPUT" ]]; then
  echo -e "\n${RED}ERROR${NC}: --out file name was not provided! Exiting..." && exit 1
fi

if [[ ! -d tmp ]]; then
  mkdir tmp 2>/dev/null
fi

echo -e "Converting transcripts IDs..."

grep '>' "$TRANSCRIPTS" | awk '{print $1,$4}' | sed 's/ gene:/_/g' > tmp/formated_IDs.lst

python3 - << END
fasta= open('$TRANSCRIPTS')
newnames= open('tmp/formated_IDs.lst')
newfasta= open('$OUTPUT', 'w')

for line in fasta:
    if line.startswith('>'):
        newname= newnames.readline()
        newfasta.write(newname)
    else:
        newfasta.write(line)

fasta.close()
newnames.close()
newfasta.close()
END

rm -R tmp/
old_ID=$(grep '>' "$TRANSCRIPTS" | head -4)
new_ID=$(grep '>' "$OUTPUT" | head -4)

echo -e "\nDONE!"
echo -e "\nThe IDs have been changed from ENSEMBL pattern:\n$old_ID\n\nTo ChimeraTE Mode 2 pattern:\n\n$new_ID"
