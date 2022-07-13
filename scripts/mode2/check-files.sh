#!/bin/bash

echo -ne "Checking files..."

if [[ -z "$MATE1" ]]; then
  echo -e "\n${RED}ERROR${NC}: --mate1 fastq files were not provided! Exiting..." && exit 0; else
    r1_samples=$(sed 's/,/\n/g' <<< "$MATE1")
    while IFS= read -r mate; do
      if [[ ! -f "$mate" ]]; then
        echo -e "\n${RED}ERROR${NC}: $mate provided with --mate1 was not found!! Exiting..." && exit 1; else
        if !(egrep -q '.fastq|.fq') <<< "$mate"; then
          echo -e "\n${RED}ERROR${NC}: Files provided with --mate1 are not in .fastq or .fq extension! Exiting..." && exit 1
        fi
      fi
    done <<< "$r1_samples"
fi

if [[ -z "$MATE2" ]]; then
  echo -e "\n${RED}ERROR${NC}: --mate2 fastq files were not provided! Exiting..." && exit 1; else
  r2_samples=$(sed 's/,/\n/g' <<< "$MATE2")
  while IFS= read -r mate; do
    if [[ ! -f "$mate" ]]; then
      echo -e "\n${RED}ERROR${NC}: $mate provided with --mate2 was not found!! Exiting..." && exit 1; else
      if !(egrep -q '.fastq|.fq') <<< "$mate"; then
        echo -e "\n${RED}ERROR${NC}: Files provided with --mate2 are not in .fastq or .fq extension! Exiting..." && exit 1
      fi
    fi
  done <<< "$r2_samples"
fi

if [[ -z "$TE_FASTA" ]]; then
  echo -e "\n${RED}ERROR${NC}: --te fasta file was not provided! Exiting..." && exit 1; else
  if !(egrep -q '.fasta|.fa') <<< "$TE_FASTA"; then
    echo -e "\n${RED}ERROR${NC}: $TE_FASTA provided with --te is not in .fasta extension! Exiting..." && exit 1; else
    if [[ ! -f "$TE_FASTA" ]]; then
      echo -e "\n${RED}ERROR${NC}: $TE_FASTA provided with --te was not found!! Exiting..." && exit 1
    fi
  fi
fi

if [[ -z "$TRANSCRIPTS_FASTA" ]]; then
  echo -e "\n${RED}ERROR${NC}: --transcripts fasta file was not provided! Exiting..." && exit 1; else
  if !(egrep -q '.fasta|.fa') <<< "$TRANSCRIPTS_FASTA"; then
    echo -e "\n${RED}ERROR${NC}: $TRANSCRIPTS_FASTA provided with --transcripts is not in .fasta extension! Exiting..." && exit 1; else
    if [[ ! -f "$TRANSCRIPTS_FASTA" ]]; then
      echo -e "\n${RED}ERROR${NC}: $TRANSCRIPTS_FASTA provided with --transcripts was not found!! Exiting..." && exit 1
    fi
  fi
fi

if [[ -z "$STRANDED" ]]; then
        echo -e "\n${RED}ERROR${NC}: The strand-specific parameter was not selected! Please use -s | --strand (rf-stranded or fr-stranded)	Exiting..." && exit 1; else
        if [ "$STRANDED" != "rf-stranded" ] && [ "$STRANDED" != "fr-stranded" ]; then
                echo -e "\n${RED}ERROR${NC}: The option -$STRANDED- is not accepted, please use --strand (rf-stranded or fr-stranded)	Exiting..." && exit 1
        fi
fi

echo -e "\t${GREEN}==> DONE!${NC}"
