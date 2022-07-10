#!/bin/bash

cut -f4 "$PROJECT"/TE-exonized_final.ct |  tr ';' '\n' | sed 's/ //g' | sort | uniq -c | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | awk '{print $2,$1}' | sed 's/ /\t/g' > "$PROJECT"/plot/TE-exon-freq.tsv

cut -f4 "$PROJECT"/TE-initiated_final.ct | grep -v '0' | sed '1d' | tr ';' '\n' | sed 's/ //g' | sort | uniq -c | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | awk '{print $2,$1}' | sed 's/ /\t/g' > "$PROJECT"/plot/TE-init-freq-UPSTREAM.tsv
cut -f8 "$PROJECT"/TE-initiated_final.ct | grep -v '0' | sed '1d' | tr ';' '\n' | sed 's/ //g' | sort | uniq -c | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | awk '{print $2,$1}' | sed 's/ /\t/g' > "$PROJECT"/plot/TE-init-freq.tsv

cut -f3 "$PROJECT"/TE-terminated_final.ct | grep -v '-' | sed '1d' | tr ';' '\n' | sed 's/ //g' | sort | uniq -c | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | awk '{print $2,$1}' | sed 's/ /\t/g' > "$PROJECT"/plot/TE-term-freq-DOWN.tsv
cut -f6 "$PROJECT"/TE-terminated_final.ct | grep -v '-' | sed '1d' | tr ';' '\n' | sed 's/ //g' | sort | uniq -c | tr -s ' ' | sed 's/^ *//g' | tr ' ' '\t' | awk '{print $2,$1}' | sed 's/ /\t/g' > "$PROJECT"/plot/TE-term-freq.tsv

cp scripts/plot.R "$PROJECT"/plot/
cd "$PROJECT"/plot/
Rscript plot.R
