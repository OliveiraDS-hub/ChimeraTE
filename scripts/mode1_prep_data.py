import argparse
import subprocess
import datetime
import glob, os
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *

if __name__ == 'mode1_prep_data':
    clock = time()
    print(str(clock) + '\t' + "Splitting genes and exons positions")
    complete_gene = pd.read_csv(str(tmp + '/gtf_file.gtf'), header=None, sep="\t", usecols=[0,1,2,3,4,5,6,7,8],names=['scaf', 'source', 'feature', 'start', 'end', 'dot', 'strand', 'dot2', 'ID'])
    complete_gene['ID'] = complete_gene['ID'].str.replace("gene-", '')
    # complete_gene['ID'] = complete_gene['ID'].str.replace("\"", '')
    # complete_gene.to_csv(str(tmp + '/gtf_file.gtf'), sep='\t', encoding='utf-8', header=None,index=False)

    # #Gene region - GTF to BED
    ID_format = complete_gene[complete_gene['ID'].str.startswith('ID=')]
    if ID_format.empty:
        complete_gene = pd.read_csv(str(tmp + '/gtf_file.gtf'), header=None, sep="\t", usecols=[0,2,3,4,5,6,7,8],names=['scaf', 'source', 'feature', 'start', 'end', 'dot', 'strand', 'dot2', 'ID'])
        complete_gene['ID'] = complete_gene['ID'].str.split(';').str[0]
        complete_gene['ID'] = complete_gene['ID'].str.replace("\"", '')
        complete_gene['ID'] = complete_gene['ID'].str.replace("gene_id ", '')
        gene = complete_gene[complete_gene["feature"] == "gene"]
        gene[['scaf', 'start', 'end', 'ID', 'dot', 'strand']].to_csv(str(tmp + '/gene_coord.bed'), sep='\t', encoding='utf-8', header=None,index=False)

        exon = complete_gene[complete_gene["feature"] == "exon"].drop_duplicates()
        exon['ID'] = exon['ID'].str.replace('.*?gene=', '').str.replace(';product=.*', '')
        exon[['scaf', 'start', 'end', 'ID', 'dot', 'strand']].to_csv(str(tmp + '/exon_file.bed'), sep='\t', encoding='utf-8', header=None,index=False)

    else:
        complete_gene = pd.read_csv(str(tmp + '/gtf_file.gtf'), header=None, sep="\t", usecols=[0,2,3,4,5,6,7,8],names=['scaf', 'source', 'feature', 'start', 'end', 'dot', 'strand', 'dot2', 'ID'])
        complete_gene['ID'] = complete_gene['ID'].str.split(';').str[0]
        complete_gene['ID'] = complete_gene['ID'].str.replace("ID=gene-", '')
        gene = complete_gene[complete_gene["feature"] == "gene"]
        gene[['scaf', 'start', 'end', 'ID', 'dot', 'strand']].to_csv(str(tmp + '/gene_coord.bed'), sep='\t', encoding='utf-8', header=None,index=False)

        complete_gene = pd.read_csv(str(tmp + '/gtf_file.gtf'), header=None, sep="\t", usecols=[0,2,3,4,5,6,7,8],names=['scaf', 'source', 'feature', 'start', 'end', 'dot', 'strand', 'dot2', 'ID'])
        exon = complete_gene[complete_gene["feature"] == "exon"].drop_duplicates()
        exon['ID'] = exon['ID'].str.replace('.*?gene=', '').str.split(';').str[0]
        exon[['scaf', 'start', 'end', 'ID', 'dot', 'strand']].drop_duplicates().to_csv(str(tmp + '/exon_file.bed'), sep='\t', encoding='utf-8', header=None,index=False)


    # TE insertions - GTF to BED
    complete_te = pd.read_csv(args.te, header=None, sep="\t", usecols=[0,2,3,4,5,6,7,8],names=['scaf', 'source', 'feature', 'start', 'end', 'dot', 'strand', 'dot2', 'ID'])
    complete_te[['scaf', 'start', 'end', 'ID', 'dot2', 'strand']].to_csv(str(tmp + '/TE_file.bed'), sep='\t', encoding='utf-8', header=None,index=False)
    print(colored("Done!", "green", attrs=['bold']))


    # Creating STAR db
    clock = time()

    if str(args.index) != "None":
        print("STAR index provided in " + str(args.index) + '\n')
    elif os.path.isfile(str(out_dir + '/index/SAindex')):
        print("Star index found! Be sure that it is not corrupted \n")
    else:
        print(str(clock) + '\t' + "Creating STAR index with " + str(out_genome))
        subprocess.call(['STAR', '--runThreadN', str(args.threads), '--runMode', str("genomeGenerate"), "--genomeDir", str(out_dir + '/index'), \
        "--genomeFastaFiles", str(args.genome), "--sjdbGTFfile", str(tmp + '/gtf_file.gtf'), "--sjdbOverhang", str(99)], stdout=subprocess.DEVNULL)
        print(colored("Done!", "green", attrs=['bold']))








#
