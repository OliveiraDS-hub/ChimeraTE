import argparse
import subprocess
import datetime
import glob, os
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *
from termcolor import colored
import datetime

def test_GTF_feature(tmp_dir, gene_GTF, feature_name):
    if not os.path.exists(f"{tmp_dir}/gtf_file.gtf") or os.stat(f"{tmp_dir}/gtf_file.gtf").st_size == 0:
        with open(f"{tmp_dir}/gtf_file.gtf", 'w') as f, open(gene_GTF, 'r') as in_gtf:
            for line in in_gtf:
                if not (line.startswith('#')):
                    f.write(line)

    gtf_file = pd.read_csv(str(f"{tmp_dir}/gtf_file.gtf"), header=None, sep="\t",usecols=[2, 8], names=['feature', 'gene_id'])

    if any(str(feature_name) in feature for feature in gtf_file['feature']):
        feature_rows = (gtf_file['feature'] == feature_name).sum()
        if feature_name == "gene":
            print(f"--> {feature_rows} {feature_name}s and", end = " ")
        else:
            print(f"{feature_rows} {feature_name}s")
        return True
    else:
        return False

def gene_IDs(tmp_dir):
    complete_gene = pd.read_csv(str(f"{tmp_dir}/gtf_file.gtf"), nrows = 1000, header=None, sep="\t", usecols=[2, 8],names=['feature', 'ID'])
    complete_gene = complete_gene[complete_gene["feature"] == "gene"]
    complete_gene['ID'] = complete_gene['ID'].str.replace("ID=gene-", '')
    complete_gene['ID'] = complete_gene['ID'].str.split(';').str[0]
    ID_format = complete_gene[complete_gene['ID'].str.startswith('ID=')]
    if ID_format.empty:
        complete_gene['ID'] = complete_gene['ID'].str.replace("\"", '').str.replace("gene_id ", '')
    else:
        complete_gene['ID'] = complete_gene['ID'].str.replace("ID=gene-", '')

    unique_genes = complete_gene['ID'].drop_duplicates()
    first_genes = unique_genes.head(2)
    first_genes = ', '.join(first_genes.tolist())
    return first_genes

def count_TE_families(TE_file):
    #from __main__ import args.te
    complete_te = pd.read_csv(str(TE_file), header=None, sep="\t", usecols=[8],names=['TE_family'])
    unique_te_families = complete_te['TE_family'].drop_duplicates()
    total_families = unique_te_families.count()
    first_TEs = unique_te_families.head(2)
    first_TEs = ', '.join(first_TEs.tolist())

    print(f"--> {total_families} TE families")
    print(f"First two TE families: {first_TEs}")
    print(colored("Done!", "green", attrs=['bold']))

def annotation_manager(tmp_dir, index, out_dir, genome, te_gtf, threads):
    clock = time()
    print(f"\n{clock}\tSplitting genes and exons positions")
    if not os.path.exists(f"{tmp_dir}/exon_file.bed"):
        complete_gene = pd.read_csv(str(f"{tmp_dir}/gtf_file.gtf"), header=None, sep="\t", usecols=[0,1,2,3,4,5,6,7,8],names=['scaf', 'source', 'feature', 'start', 'end', 'dot', 'strand', 'dot2', 'ID'])
        complete_gene['ID'] = complete_gene['ID'].str.replace("ID=gene-", '')

        ### Gene region - GTF to BED
        ########complete_gene = pd.read_csv(str(f"{tmp}/gtf_file.gtf"), header=None, sep="\t", usecols=[0,2,3,4,5,6,7,8],names=['scaf', 'source', 'feature', 'start', 'end', 'dot', 'strand', 'dot2', 'ID'])
        complete_gene['ID'] = complete_gene['ID'].str.split(';').str[0]
        ID_format = complete_gene[complete_gene['ID'].str.startswith('ID=')]

        if ID_format.empty:
            complete_gene['ID'] = (
                complete_gene['ID']
                .str.replace('"', '', regex=False)
                .str.replace('gene_id ', '', regex=False)
                .str.replace('transcript_id ', '', regex=False)
                .str.replace(r'\.t\d+', '', regex=True)
            )            
            gene = complete_gene[complete_gene["feature"] == "gene"]
            gene[['scaf', 'start', 'end', 'ID', 'dot', 'strand']].to_csv(f"{tmp_dir}/gene_coord.bed", sep='\t', encoding='utf-8', header=None,index=False)

            exon = complete_gene[complete_gene["feature"] == "exon"].drop_duplicates()
            exon['ID'] = exon['ID'].str.replace('.*?gene=', '').str.replace(';product=.*', '')
            exon[['scaf', 'start', 'end', 'ID', 'dot', 'strand']].to_csv(f"{tmp_dir}/exon_file.bed", sep='\t', encoding='utf-8', header=None,index=False)
        else:
            complete_gene['ID'] = complete_gene['ID'].str.replace("ID=gene-", '')
            gene = complete_gene[complete_gene["feature"] == "gene"]
            gene[['scaf', 'start', 'end', 'ID', 'dot', 'strand']].to_csv(f"{tmp_dir}/gene_coord.bed", sep='\t', encoding='utf-8', header=None,index=False)

            complete_gene = pd.read_csv(str(tmp_dir + '/gtf_file.gtf'), header=None, sep="\t", usecols=[0,2,3,4,5,6,7,8],names=['scaf', 'source', 'feature', 'start', 'end', 'dot', 'strand', 'dot2', 'ID'])
            exon = complete_gene[complete_gene["feature"] == "exon"].drop_duplicates()
            exon['ID'] = exon['ID'].str.replace('.*?gene=', '').str.split(';').str[0]
            exon[['scaf', 'start', 'end', 'ID', 'dot', 'strand']].drop_duplicates().to_csv(f"{tmp_dir}/exon_file.bed", sep='\t', encoding='utf-8', header=None,index=False)


    ### TE insertions - GTF to BED
    if not os.path.exists(f'{tmp_dir}/TE_file.bed'):
        complete_te = pd.read_csv(te_gtf, header=None, sep="\t", usecols=[0,2,3,4,5,6,7,8],names=['scaf', 'source', 'feature', 'start', 'end', 'dot', 'strand', 'dot2', 'ID'])
        complete_te[['scaf', 'start', 'end', 'ID', 'dot2', 'strand']].to_csv(f'{tmp_dir}/TE_file.bed', sep='\t', encoding='utf-8', header=None,index=False)
        print(colored("Done!", "green", attrs=['bold']))


    ### Creating STAR db
    clock = time()
    if str(index) != "None":
        print(f"STAR index provided in {index}\n")
    elif os.path.isfile(str(f"{out_dir}/index/SAindex")):
        print("STAR index found! Make sure it is not corrupted \n")
    else:
        print(f"{clock}\tCreating STAR index with {genome}")
        subprocess.call(['STAR', '--runThreadN', str(threads), '--runMode', str("genomeGenerate"), "--genomeDir", str(f"{out_dir}/index"), \
        "--genomeFastaFiles", str(genome), "--sjdbGTFfile", str(f"{tmp_dir}/gtf_file.gtf"), "--sjdbOverhang", str(99)], stdout=subprocess.DEVNULL)
        print(colored("Done!", "green", attrs=['bold']))








#
