#!/bin/bash
set -e

function usage() {
cat << help

Conversion of NCBI native transcript IDs to the ChimeraTE Mode 2 format

#Mandatory arguments:

  --transcripts     transcripts downloaded from NCBI (.fasta)
                    i.e.: GCF_000001735.4_TAIR10.1_rna.fna  (A. thaliana)
                          GCF_000001215.4_Release_6_plus_ISO1_MT_rna.fna (D. melanogaster)

  --out             output file name (.fasta)
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

grep '>' "$TRANSCRIPTS" | awk '{print $1}' | sed 's/_/-/g' > tmp/part1.rm
grep '>' "$TRANSCRIPTS" | sed 's/>.*(//g; s/).*//' | sed 's/_/-/g' | sed 's/ /-/g' > tmp/part2.rm

paste tmp/part1.rm tmp/part2.rm --delimiter '_' > tmp/formated_IDs.lst

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
echo -e "\nThe IDs have been changed from NCBI pattern:\n$old_ID\n\nTo ChimeraTE Mode 2 pattern:\n\n$new_ID"
