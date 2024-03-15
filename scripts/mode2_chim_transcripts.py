from __main__ import *
import argparse
from multiprocessing.dummy import Pool
import datetime
import subprocess
import glob
import pandas as pd
import re
import csv
from io import StringIO
import pybedtools
import contextlib
import sys
import os
import math

def prep_data():
    from __main__ import aln_dir
    pd.read_table(str(f"{aln_dir}/genes_chim.bed"), sep="\t", usecols=[0], names=['id']).drop_duplicates().to_csv(str(f"{aln_dir}/gene.lst"), header=None, index=False, sep="\t")

def chim_transcript(gene_list):
    from __main__ import aln_dir
    from __main__ import out_group

    check = None
    chim_genes_complete = pd.read_csv(str(f"{aln_dir}/genes_chim.bed"), sep="\t", usecols=[0,1,2,3,4,5], names=['transc_id', 'transc_start', 'transc_end', 'transc_read', 'transc_dot', 'transc_strand'])

    TE_chim = pd.read_csv(str(f"{aln_dir}/TEs_chim.bed"), sep="\t", usecols=[0,1,2,3,4,5], names=['TE_id', 'TE_start', 'TE_end', 'TE_read', 'TE_dot', 'TE_strand'])

    gene_id = gene_list
    gene_data = chim_genes_complete.query('transc_id == @gene_id')
    read_gene_ids = gene_data.iloc[:, 3].drop_duplicates().to_csv(header=None, index=False, sep="\t").splitlines()

    check = TE_chim[TE_chim['TE_read'].isin(read_gene_ids)].iloc[:, 0].value_counts().rename_axis('TE_fam').reset_index(name='counts').sort_values(by=['counts'], ascending=False).head(1)
    if check is not None:
        TE_family = check.loc[:, 'TE_fam'].to_csv(header=None, index=False, sep=" ").strip('\n')
        cov = check.loc[:, 'counts'].replace('\n', '', regex=True).to_csv(header=None, index=False, sep=" ").strip('\n')
        if cov and int(cov) > 0:
            with open(str(f"{out_group}/chimTEs_raw.tsv"), "a") as f:
                print(str(f"{gene_id}\t{TE_family}\t{cov}"), file=f)

def merging_transc():
    from __main__ import aln_dir
    from __main__ import out_group

    clock = time()
    print(str(f"{clock}\tMerging coverage from different isoforms..."))
    if check_file(f"{out_group}/chimTEs_final.tsv") == False:
        chimeras = pd.read_csv(str(f"{out_group}/chimTEs_raw.tsv"), sep="\t", usecols=[0,1,2], names=['gene_transc', 'TE_id', 'cov'])
        chimeras[['isoform','gene']] = chimeras.gene_transc.str.split("_",expand=True)
        gene_ids = chimeras.loc[:, 'gene'].drop_duplicates().to_list()

        for gene_id in gene_ids:
            chim_data = chimeras.query('gene == @gene_id')
            TE_fam = chim_data['TE_id'].head(1).to_csv(header=None, index=False, sep=" ").strip('\n')
            total_cov = chim_data['cov'].sum()
            isoforms = chim_data['isoform'].to_csv(header=None, index=False, sep=" ").replace('\n', str('_' + gene_id + '; '))
            isoforms = isoforms[:-2]
            with open(str(f"{out_group}/chimTEs_final.tsv"), "a") as f:
                print(str(gene_id + '\t' + TE_fam + '\t' + str(total_cov) + '\t' + str(isoforms)), file = f)
            f.close()
        print(colored("Done!", "green", attrs=['bold']))
    else:
        print(f"Merged file has been found!"); print(colored("Skipping...\n", "yellow", attrs=['bold']))

def expression():
    from __main__ import aln_dir
    from __main__ import out_group
    from __main__ import group
    from __main__ import tmp

    chim_complete = pd.read_csv(str(f"{out_group}/chimTEs_final.tsv"), sep="\t", usecols=[0,1,2,3], names=['gene_id', 'TE_id', 'cov', 'isoforms'])
    chim_gene_id = pd.read_csv(str(f"{out_group}/chimTEs_final.tsv"), sep="\t", usecols=[0], names=['col1'])
    chim_gene_id = chim_gene_id['col1'].tolist()
    transcripts_list = pd.read_csv(str(f"{out_group}/chimTEs_final.tsv"), sep="\t", usecols=[3], names=['transc_id'])
    transcripts_list = transcripts_list['transc_id'].to_csv(header=None, index=False, sep=" ").replace('; ', '\n').replace('"', '')
    transcripts_list = StringIO(transcripts_list)
    transcripts_df = pd.read_table(transcripts_list, sep="\t", usecols=[0], names=['transc_gene'])
    transcripts_df[['isoform','gene']] = transcripts_df.transc_gene.str.split("_",expand=True)

    if check_file(str(f"{aln_dir}/fpkm_counts/results.xprs")) == True:
        fpkm_gene = pd.read_csv(str(f"{aln_dir}/fpkm_counts/results.xprs"), sep="\t", usecols=[1,12], names=['isoform_id', 'fpkm'])
    else:
        fpkm_gene = pd.read_csv(str(f"{aln_dir}/genes_total_expressed.bed"), sep="\t", usecols=[0], names=['isoform_id'])

    ncol_fpkm_gene = len(fpkm_gene.columns)
    for row in chim_gene_id:
        exp_level = ''
        gene_data = transcripts_df.query('gene == @row')

        all_info = chim_complete.query('gene_id == @row')
        TE_fam = all_info.iloc[:, 1].to_csv(header=None, index=False, sep=" ").strip('\n')
        cov = all_info.iloc[:, 2].to_csv(header=None, index=False, sep=" ").strip('\n')
        isoforms = all_info.iloc[:, 3].to_csv(header=None, index=False, sep=" ").strip('\n').replace('"', '')


        for transc in gene_data.itertuples():
            transcript = transc[1]
            if ncol_fpkm_gene > 1:
                exp_level += fpkm_gene.query('isoform_id == @transcript').loc[:, ['fpkm']].to_csv(header=None, index=False, sep=" ")

        if ncol_fpkm_gene > 1:
            freq = None
            freq = exp_level.count('\n')

            if freq > 1:
                exp_level = exp_level.splitlines()
                floating = [float(x) for x in exp_level]
                total = math.fsum(floating)
                avg = total / freq
            else:
                avg = str(exp_level).strip('\n')
        else:
            avg = "NA"

        with open(str(out_group + '/' + group + '_chimreads_evidence.tsv'), "a") as f:
            print(str(row + '\t' + TE_fam + '\t' + str(cov) + '\t' + str(isoforms) + '\t' + str(avg)), file = f)
        f.close()

    copy(str(out_group + '/' + group + '_chimreads_evidence.tsv'), tmp)

def multicore_chimeras():
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import out_group

    clock = time()
    print(str(clock) + '\t' + "Identifying chimeric transcripts...")

    with open(str(aln_dir + '/gene.lst')) as f:
        gene_list = f.read().splitlines()
    f.close
    if check_file(f"{out_group}/chimTEs_raw.tsv") == False:
        pool = Pool(processes=int(args.threads))
        pool.map(chim_transcript, gene_list)
        pool.close()
        print(colored("Done!", "green", attrs=['bold']))
    else:
        print(f"Chimeric transcripts file has been found!"); print(colored("Skipping...\n", "yellow", attrs=['bold']))
