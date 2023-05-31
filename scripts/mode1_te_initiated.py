import multiprocessing
import concurrent.futures
import argparse
import datetime
import subprocess
from multiprocessing.dummy import Pool
from termcolor import colored
import glob, os
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *

def te_init(aln_dir, group, out_group):
    print("\n###########################\n## TE-initiated analysis ##\n###########################\n")
    te_init=open(str(out_group + "/TE-initiated-" + str(group) + ".tsv"), 'w')
    clock = time()
    print(str(clock) + '\t' + "Generating upstream region for expressed genes...")
    up_window = pd.read_table(str(aln_dir + '/genes_total_expressed.bed'), header=None, usecols=[0,1,2,3,4,5],names=['gene_id', 'start', 'end', 'id', 'dot', 'strand'])

    with open(str(aln_dir + '/genes_UP_window.bed'), encoding='utf-8', mode='a') as up_window_file:
        for index, row in up_window.iterrows():
            scaf = repr(row['gene_id'])
            start = repr(row['start'])
            end = repr(row['end'])
            part = repr(row['id'])
            dot = repr(row['dot'])
            strand = row['strand']
            if strand == "+":
                num = (int(start) - int(args.window))
                if num < 1:
                    num = 1
                up_window_file.write(scaf.replace("'","") + "\t" + str(num) + "\t" + str(start) + "\t" + part.replace("'","") + "\t" + dot.replace("'","") + "\t" + strand.replace("'","") + "\n")
            else:
                num = (int(end) + int(args.window))
                up_window_file.write(scaf.replace("'", "") + "\t" + str(end) + "\t" + str(num) + "\t" + part.replace("'","") + "\t" + dot.replace("'", "") + "\t" + strand.replace("'", "") + "\n")
    up_window_file.close()

    genes_UP_window = pybedtools.BedTool(str(aln_dir + '/genes_UP_window.bed'))
    chimeric_TEs = pybedtools.BedTool(str(aln_dir + '/chimeric_TEs.bed'))
    genes_TE_UP_mbed = chimeric_TEs.intersect(genes_UP_window, wa=True, wb=True, nonamecheck=True)
    bed2string = str(genes_TE_UP_mbed)
    bed2string = StringIO(bed2string)
    global df
    df = pd.read_csv(bed2string, sep="\t", header=None, usecols=[0,1,2,3,4,5,6,7,8,9,10,11],names=['scaf', 'start', 'end', 'ID', 'dot', 'strand','scaf2', 'start2', 'end2', 'ID2', 'dot2', 'strand2']).drop_duplicates()
    df["ID2"].drop_duplicates().to_csv(str(aln_dir + '/init_chimeras.lst'), encoding='utf-8', header=None,index=False)

    print(colored("Done!", "green", attrs=['bold']))

def init_chimeras(init_list):
    from __main__ import aln_dir
    gene_reads_fwd = pybedtools.BedTool(aln_dir + '/gene_reads_fwd.bed')
    TE_reads_fwd = pybedtools.BedTool(aln_dir + '/TE_reads_fwd.bed')
    gene_reads_rev = pybedtools.BedTool(aln_dir + '/gene_reads_rev.bed')
    TE_reads_rev = pybedtools.BedTool(aln_dir + '/TE_reads_rev.bed')

    exon = pd.read_csv(str(aln_dir + '/chim_exons.bed'), sep="\t", usecols=[0,1,2,3,4,5], names=['scaf_exon', 'start_exon', 'end_exon', 'ID_exon', 'dot_exon', 'strand_exon'])

    gene_id = init_list
    TEs_upstream = df.query('ID2 == @gene_id')
    gene_coord = df.query('ID2 == @gene_id').iloc[:, 6:12].drop_duplicates()
    chr_gene = gene_coord['scaf2'].to_string(index=False).replace(" ","")
    s_gene = gene_coord['start2'].to_string(index=False).replace(" ","")
    e_gene = gene_coord['end2'].to_string(index=False).replace(" ","")
    strand = gene_coord['strand2'].to_string(index=False).replace(" ","")
    exon_coord = exon.query('ID_exon == @gene_id').to_csv(sep='\t', encoding='utf-8', header=None,index=False)
    exon_bed = pybedtools.BedTool(str(exon_coord), from_string=True)

    if exon_coord:
        exons_coord = StringIO(exon_coord)
        strand_exon = pd.read_csv(exons_coord, sep = "\t", header=None, usecols=[5]).iloc[0].item()

        for row in TEs_upstream.itertuples():
            TE_tmp = row[1:7]
            TE_tmp = pd.DataFrame(TE_tmp).transpose()

            chr_TE = row[1]
            s_TE = row[2]
            e_TE = row[3]
            TE_family = row[4]
            dot = row[5]
            TE_strd = row[6]

            merging_TE_reads = None
            merging_gene_reads = None

            if strand == "+":
                if int(e_TE) > int(e_gene):
                    e_TE = int(e_gene)
                TE_coord = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
                TE_bed = pybedtools.BedTool(TE_coord, from_string=True)
                read_genes = gene_reads_fwd.intersect(exon_bed, wa=True, wb=True, nonamecheck=True)
                intersect_gene = str(read_genes)
                intersect_gene = StringIO(intersect_gene)
                reads_gene_col = pd.read_table(intersect_gene, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

                reads_TE = TE_reads_fwd.intersect(TE_bed, wa=True, wb=True, nonamecheck=True, f=str(args.overlap))
                intersect_TE = str(reads_TE)
                intersect_TE = StringIO(intersect_TE)
                reads_TE_col = pd.read_table(intersect_TE, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()
            else:
                if int(s_TE) < int(s_gene):
                        s_TE = int(s_gene)
                TE_coord = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
                TE_bed = pybedtools.BedTool(TE_coord, from_string=True)
                read_genes = gene_reads_rev.intersect(exon_bed, wa=True, wb=True, nonamecheck=True)
                intersect_gene = str(read_genes)
                intersect_gene = StringIO(intersect_gene)
                reads_gene_col = pd.read_table(intersect_gene, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

                reads_TE = TE_reads_rev.intersect(TE_bed, wa=True, wb=True, nonamecheck=True, f=str(args.overlap))
                intersect_TE = str(reads_TE)
                intersect_TE = StringIO(intersect_TE)
                reads_TE_col = pd.read_table(intersect_TE, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

            if not reads_gene_col.empty:
                if merging_gene_reads is not None:
                    merging_gene_reads = merging_gene_reads.append(reads_gene_col, ignore_index=True).drop_duplicates()
                else:
                    merging_gene_reads = reads_gene_col

            if not reads_TE_col.empty:
                if merging_TE_reads is not None:
                    merging_TE_reads = merging_TE_reads.append(reads_TE_col, ignore_index=True).drop_duplicates()
                else:
                    merging_TE_reads = reads_TE_col

            if merging_TE_reads is not None:
                if merging_gene_reads is not None:
                    list_TE_reads = merging_TE_reads['read_ID'].tolist()
                    match = merging_gene_reads[merging_gene_reads['read_ID'].isin(list_TE_reads)].drop_duplicates()
                    cov = len(match)
                    if cov > 0:
                        with open(str("TE-initiated_tmp.tsv"), "a") as te_init:
                            print(str(gene_id) + "\t" + str(strand_exon) + "\t" + str(chr_gene) + ":" + str(s_gene) + "-" + str(e_gene) + "\t" + str(TE_family) + "\t" + str(TE_strd) + "\t" + str(chr_TE) + ":" + str(s_TE) + "-" + str(e_TE) + "\t" + str(cov), file=te_init)
                        te_init.close()

def multicore_process_init():
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import out_group

    clock = time()
    print(str(clock) + '\t' + "Searching for TE-initiated transcripts...")

    #TE-initiated chimeras
    with open(str(aln_dir + '/init_chimeras.lst')) as f:
        init_list = f.read().splitlines()
    f.close
    pool = Pool(processes=int(args.threads))
    pool.map(init_chimeras, init_list)
    pool.close()

    init_output = str('TE-initiated_tmp.tsv')
    if os.path.exists(init_output) == True:
        if os.stat(str('TE-initiated_tmp.tsv')).st_size > 0:
            os.rename('TE-initiated_tmp.tsv', str('TE-initiated-' + str(group) + '.tsv'))
            copy(str('TE-initiated-' + str(group) + '.tsv'), out_group)
            os.remove('TE-initiated-' + str(group) + '.tsv')
            print(colored("Done!", "green", attrs=['bold']))
            pybedtools.cleanup(remove_all=True)
    else:
        print(colored('There are no TE-initiated transcripts!', "red"))
