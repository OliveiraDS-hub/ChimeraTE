import argparse
import datetime
import subprocess
import glob, os
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *

def replicability(chim_type):
    df_replicates = pd.DataFrame()
    clock = time()
    print(f"{clock}\tSearching for {chim_type} transcripts found in at least {args.replicate} replicates")
    for file in os.listdir(str(tmp)):
        if file.startswith(str(chim_type)):
            if os.stat(str(f"{tmp}/{file}")).st_size > 0:
                data = pd.read_csv(str(f"{tmp}/{file}"), sep = '\t', header = None)
                df_replicates = df_replicates.append(data, ignore_index=True)
    return df_replicates

list = ["TE-initiated", "TE-terminated", "TE-exonized"]

for chimera in list:
    replicated = replicability(chimera)
    if not replicated.empty:
        if chimera == "TE-exonized":
            duplicate_mask = replicated.duplicated(subset=[0, 3, 5, 7], keep=False)
            # Find rows that are duplicated at least minimum of --replicate based on columns 0 (gene_id), 3 (TE_id), and 5 (TE_position)
            filtered_df = replicated[duplicate_mask].groupby([0, 3, 5, 7]).filter(lambda x: len(x) >= int(args.replicate))
            if not filtered_df.empty:
                ### Calculate the mean of chimeric reads on chimeras found in multiple replicates
                mean_values = filtered_df.groupby([0,1,2,3,4,5,7])[6].mean().round(2).reset_index()

                ### Select chimeras with mean greater than minimum coverage
                mean_values = mean_values[mean_values[6] >= int(args.coverage)]
                mean_values.columns = ['gene_id', 'gene_strand', 'gene_pos', 'TE_id', 'TE_strand', 'TE_pos', 'exonized_type', 'chim_reads']
                mean_values.to_csv(f"{out_dir}/{chimera}_final.tsv", sep='\t', header= True, index=False)
                row_count = len(mean_values)
                print(f"{row_count} {chimera} transcripts have been detected!")
            else:
                print(f"There is no {chimera} transcript replicated at least {args.replicate} times\n")
        else:
            duplicate_mask = replicated.duplicated(subset=[0, 3, 5], keep=False)
            # Find rows that are duplicated at least minimum of --replicate based on columns 0 (gene_id), 3 (TE_id), and 5 (TE_position)
            filtered_df = replicated[duplicate_mask].groupby([0, 3, 5]).filter(lambda x: len(x) >= int(args.replicate))
            if not filtered_df.empty:
                mean_values = filtered_df.groupby([0, 1, 2, 3, 4, 5])[6].mean().round(2).reset_index()
                mean_values = mean_values[mean_values[6] >= int(args.coverage)]
                mean_values.columns = ['gene_id', 'gene_strand', 'gene_pos', 'TE_id', 'TE_strand', 'TE_pos', 'chim_reads']
                mean_values.to_csv(f"{out_dir}/{chimera}_final.tsv", sep='\t', header=True, index=False)
                row_count = len(mean_values)
                print(f"{row_count} {chimera} transcripts have been detected!")
            else:
                print(f"There is no {chimera} transcript replicated at least {args.replicate} times\n")
        print(colored("Done!", "green", attrs=['bold']))
    else:
        print(f"There is no {chimera} transcripts detected.\n")








#