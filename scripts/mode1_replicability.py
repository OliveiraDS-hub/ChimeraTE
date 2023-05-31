import argparse
import datetime
import subprocess
import glob, os
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *

replicate = str(args.replicate)
list = ["TE-initiated", "TE-terminated", "TE-exonized"]

for chim_type in list:
    data_you_need=pd.DataFrame()
    clock = time()
    print(str(clock) + '\t' + "Searching for " + str(chim_type) + " transcripts found in at least " + str(replicate) + " replicates" )
    for file in os.listdir(str(tmp)):
        if file.startswith(str(chim_type)):
            table = str(tmp + '/' + file)
            if os.path.exists(table) == True:
                if os.stat(str(table)).st_size > 0:
                    data = pd.read_csv(str(tmp + "/" + file), sep = '\t', header = None)
                    data_you_need=data_you_need.append(data, ignore_index=True)

    if data_you_need is not None and not data_you_need.empty:
        data_you_need.iloc[:,[0,3,5]].replace('"', '').to_csv(str(tmp) +  '/' + "merging.tsv", sep='&', header=None, index=False)
        pd.set_option('display.max_colwidth', None)
        merged_table = pd.read_csv(str(tmp) +  '/' + "merging.tsv", sep = "\t", usecols=[0],names=['chimeras']).replace('"', '', regex=True).value_counts().rename_axis('chimeras').reset_index(name='counts')
        merged_table = merged_table.query('counts == @replicate').iloc[:,0]#.to_string(header=False, index=False).replace('&', '\t').replace(' ', '')

        if not merged_table.empty:
            merged_table = merged_table.to_string(header=False, index=False).replace('&', '\t').replace(' ', '')
            merged_table = pd.DataFrame([x.split('\t') for x in merged_table.split('\n')])
            output = None
            for row in merged_table.itertuples():
                merged = None
                gene_id = row[1]
                TE_family = row[2]
                TE_pos = row[3]

                for file in os.listdir(str(tmp)):
                    if file.startswith(str(chim_type)):
                        if chim_type == "TE-exonized":
                            df = pd.read_csv(str(tmp + "/" + file), sep = '\t', usecols=[0,1,2,3,4,5,6,7],names=['gene_id', 'gene_strand', 'gene_pos', 'TE_id', 'TE_strand', 'TE_pos', 'chim_reads', 'exonized_type'])
                        else:
                            df = pd.read_csv(str(tmp + "/" + file), sep = '\t', usecols=[0,1,2,3,4,5,6],names=['gene_id', 'gene_strand', 'gene_pos', 'TE_id', 'TE_strand', 'TE_pos', 'chim_reads'])

                        match = df.query('gene_id == @gene_id').query('TE_id == @TE_family').query('TE_pos == @TE_pos')
                        if merged is None:
                            merged = match
                        else:
                            merged = merged.append(match, ignore_index=True)
                if chim_type == "TE-exonized":
                    chim_info = merged[['gene_id', 'gene_strand', 'gene_pos', 'TE_id', 'TE_strand', 'TE_pos', 'exonized_type']].drop_duplicates()
                else:
                    chim_info = merged[['gene_id', 'gene_strand', 'gene_pos', 'TE_id', 'TE_strand', 'TE_pos']].drop_duplicates()

                mean = merged["chim_reads"].mean()
                if mean >= int(args.coverage):
                    chim_info['chim_reads'] = mean
                    if output is None:
                        output = chim_info
                    else:
                        output = output.append(chim_info, ignore_index=True)
            if output is None:
                print("There is no " + str(chim_type) + " transcripts found in at least " + str(replicate) + " replicates!")
            else:
                output.to_csv(str(out_dir) + '/' + str(chim_type) + '_final.tsv', sep='\t', header= True, index=False)
            os.remove(str(tmp) +  '/' + "merging.tsv")
        else:
            print("There is no " + str(chim_type) + " transcripts found in at least " + str(replicate) + " replicates!")
    else:
        print("There is no " + str(chim_type) + " transcripts found in at least " + str(replicate) + " replicates!")
    print(colored("Done!", "green", attrs=['bold']))
