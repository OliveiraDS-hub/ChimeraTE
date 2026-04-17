import argparse
import datetime
import subprocess
import glob, os
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *
from pathlib import Path

def merge_dfs(chim_type, rep_n):
    from __main__ import tmp
    df_replicates = pd.DataFrame()
    clock = time()
    print(f"{clock}\tSearching for {chim_type} transcripts found in at least {rep_n} replicates")

    files = list(Path(tmp).glob(f"{chim_type}*.tsv"))
    out_dfs = [
        pd.read_csv(f, sep='\t', header=None)
        for f in files
        if f.stat().st_size > 0]

    if out_dfs:  # only concat if there is at least one DataFrame
        df_replicates = pd.concat(out_dfs, ignore_index=True)
        return df_replicates
    else:
        print(f"No {chim_type} results for at least {rep_n}!")
        return None
        

def replicability(rep_n, coverage):
    from __main__ import out_dir

    for chimera in ["TE-initiated", "TE-terminated", "TE-exonized_embedded", "TE-exonized_overlapped", "TE-exonized_intronic"]:
        replicated = merge_dfs(chimera, rep_n)
        if replicated is not None:
            # if chimera == "TE-exonized":
            #     duplicate_mask = replicated.duplicated(subset=[0, 3, 5, 7], keep=False)
            #     # Find rows that are duplicated at least minimum of --replicate based on columns 0 (gene_id), 3 (TE_id), and 5 (TE_position)
            #     filtered_df = replicated[duplicate_mask].groupby([0, 3, 5, 7]).filter(lambda x: len(x) >= int(rep_n))
            #     if not filtered_df.empty:
            #         ### Calculate the mean of chimeric reads on chimeras found in multiple replicates
            #         mean_values = filtered_df.groupby([0,1,2,3,4,5,7])[6].mean().round(2).reset_index()

            #         ### Select chimeras with mean greater than minimum coverage
            #         mean_values = mean_values[mean_values[6] >= int(coverage)]
            #         mean_values.columns = ['gene_id', 'gene_strand', 'gene_pos', 'TE_id', 'TE_strand', 'TE_pos', 'exonized_type', 'chim_reads']
            #         mean_values.to_csv(f"{out_dir}/{chimera}_final.tsv", sep='\t', header= True, index=False)
            #         row_count = len(mean_values)
            #         print(f"{row_count} {chimera} transcripts have been detected!")
            #     else:
            #         print(f"There is no {chimera} transcript replicated at least {rep_n} times\n")
            # else:
            duplicate_mask = replicated.duplicated(subset=[0, 3, 5], keep=False)
            # Find rows that are duplicated at least minimum of --replicate based on columns 0 (gene_id), 3 (TE_id), and 5 (TE_position)
            filtered_df = replicated[duplicate_mask].groupby([0, 3, 5]).filter(lambda x: len(x) >= int(rep_n))
            if not filtered_df.empty:
                # print(filtered_df)
                mean_values = filtered_df.groupby([0, 1, 2, 3, 4, 5])[6].mean().round(2).reset_index()
                mean_values = mean_values[mean_values[6] >= int(coverage)]
                mean_values.columns = ['gene_id', 'gene_strand', 'gene_pos', 'TE_id', 'TE_strand', 'TE_pos', 'chim_reads']
                mean_values.to_csv(f"{out_dir}/{chimera}_final.tsv", sep='\t', header=True, index=False)
                row_count = len(mean_values)
                print(f"{row_count} {chimera} transcripts have been detected!")
            else:
                print(f"There is no {chimera} transcript replicated at least {rep_n} times\n")
            print(colored("Done!", "green", attrs=['bold']))








