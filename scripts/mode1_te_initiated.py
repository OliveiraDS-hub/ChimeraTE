from multiprocessing import Pool
from termcolor import colored
import concurrent.futures
import argparse
import datetime
import subprocess
import glob, os
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *
import shutil

def te_init(aln_dir, window):
    print("\n###########################\n## TE-initiated analysis ##\n###########################\n")

    clock = time()
    print(f'{clock}\tGenerating upstream region for expressed genes...')
    def create_promoter_bed(input_bed, output_bed, window):
        ### Function to create promoter region based on the windows size (3k default) | output = genes_UP_window.bed
        df = pd.read_csv(input_bed, sep="\t", header=None, names=['gene_id', 'start', 'end', 'id', 'dot', 'strand'])

        # Create new columns for promoter start and end
        df['prom_start'] = df.apply(
            lambda row: max(1, row['start'] - window) if row['strand'] == '+' else row['end'],
            axis=1
        )
        df['prom_end'] = df.apply(
            lambda row: row['start'] if row['strand'] == '+' else row['end'] + window,
            axis=1
        )

        # Reorder columns for BED format
        result = df[['gene_id', 'prom_start', 'prom_end', 'id', 'dot', 'strand']]

        # Save as BED
        result.to_csv(output_bed, sep="\t", header=False, index=False)
    ### Create bed file with promoter regions | output = genes_UP_window.bed
    create_promoter_bed(f'{aln_dir}/genes_total_expressed.bed', f'{aln_dir}/genes_UP_window.bed', window)
 
    ### Intersect promoters with TEs and get the list of the genes | output = init_chimeras.lst
    genes_TE_UP = intersection_any_bp_mbed(f'{aln_dir}/genes_UP_window.bed', f'{aln_dir}/chimeric_TEs.bed') 
    genes_TE_UP_IDs = get_IDs_from_bed(genes_TE_UP)
    genes_TE_UP_IDs["Name"].drop_duplicates().to_csv(f'{aln_dir}/init_chimeras.lst', header=None, index=False)
    
    print(colored("Done!", "green", attrs=['bold']))
    return genes_TE_UP

def init_chimeras(gene_id, itsct_mbed_TEs_up):
    from __main__ import aln_dir

    ### Extract all rows (TEs itersecting promoters)
    TEs_upstream = itsct_mbed_TEs_up.query('Name == @gene_id')[['Chromosome', 'Start_b', 'End_b', 'Name_b', 'Score_b', 'Strand_b']].drop_duplicates()
    TEs_upstream.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    
    ### Extract the promoter region - 6 columns
    gene_coord = itsct_mbed_TEs_up.query('Name == @gene_id').iloc[:, 0:6].drop_duplicates()
    ### Extract only chr
    chr_gene = gene_coord['Chromosome'].values[0]
    ### Extract only start
    s_gene   = gene_coord['Start'].values[0]
    ### Extract only end
    e_gene   = gene_coord['End'].values[0]
    ### Extract only strand
    strand   = gene_coord['Strand'].values[0]



    ### Load chimeric exons
    exon = pd.read_csv(f'{aln_dir}/chim_exons.bed', sep="\t", names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'])
    ### Extract all exons for this gene
    exon_coord = exon.query('Name == @gene_id')#.to_csv(sep='\t', encoding='utf-8', header=None,index=False)

    if not exon_coord.empty:        
        for row in TEs_upstream.itertuples():
            TE_coord_bed = pd.DataFrame([row._asdict()]).drop(columns='Index', errors='ignore')
            s_TE = row[2]
            e_TE = row[3]
            
            chim_reads_TE = pd.DataFrame()
            chim_reads_exons = pd.DataFrame()

            if strand == "+":
                if int(e_TE) > int(e_gene):
                    TE_coord_bed.at[TE_coord_bed.index[0], "End"] = int(e_gene)
                ### Check if TE located upstream has chimeric reads
                chim_reads_TE = intersection_any_bp(f'{aln_dir}/TE_reads_fwd.bed', TE_coord_bed)
                if isinstance(chim_reads_TE, pd.DataFrame):
                    ### Extract chimeric reads from TE
                    chim_reads_TE = get_IDs_from_bed(chim_reads_TE)
                    chim_reads_TE = chim_reads_TE[['Name']].replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

                    ### Extract chimeric reads from exon
                    chim_reads_exons = intersection_any_bp(f'{aln_dir}/gene_reads_rev.bed', exon_coord) 
                    if chim_reads_exons is not None:
                        chim_reads_exons = get_IDs_from_bed(chim_reads_exons)
                        chim_reads_exons = chim_reads_exons[['Name']].replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()
            
            else:
                if int(s_TE) < int(s_gene):
                    TE_coord_bed.at[TE_coord_bed.index[0], "Start"] = int(s_gene)
                ### Check if TE located upstream has chimeric reads
                chim_reads_TE = intersection_any_bp(f'{aln_dir}/TE_reads_rev.bed', TE_coord_bed)
                if isinstance(chim_reads_TE, pd.DataFrame):
                    ### Extract chimeric reads from TE
                    chim_reads_TE = get_IDs_from_bed(chim_reads_TE)
                    chim_reads_TE = chim_reads_TE[['Name']].replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

                    ### Extract chimeric reads from exon
                    chim_reads_exons = intersection_any_bp(f'{aln_dir}/gene_reads_rev.bed', exon_coord) 
                    if chim_reads_exons is not None:
                        chim_reads_exons = get_IDs_from_bed(chim_reads_exons)
                        chim_reads_exons = chim_reads_exons[['Name']].replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()
            
            if (chim_reads_TE is not None and not chim_reads_TE.empty) and \
                (chim_reads_exons is not None and not chim_reads_exons.empty):
                cov = chim_reads_TE["Name"].isin(chim_reads_exons["Name"]).sum()
                if cov > 0:
                    chr_TE = row[1]
                    TE_strd = row[6]
                    TE_family = row[4]
                    strand_exon = str(exon_coord.iloc[0]["Strand"])
                    # print(f'{gene_id}\t{strand_exon}\t{chr_gene}:{s_gene}-{e_gene}\t{TE_family}\t{TE_strd}\t{chr_TE}:{s_TE}-{e_TE}\t{cov}')
                    with open(f'{aln_dir}/TE-initiated_tmp.tsv', 'a') as te_init:
                        te_init.write(f'{gene_id}\t{strand_exon}\t{chr_gene}:{s_gene}-{e_gene}\t{TE_family}\t{TE_strd}\t{chr_TE}:{s_TE}-{e_TE}\t{cov}\n')
                    te_init.close()

def multicore_process_init(itsct_mbed_TEs_up, threads):
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import tmp

    clock = time()
    print(f'{clock}\tSearching for TE-initiated transcripts...')

    #TE-initiated chimeras
    with open(f'{aln_dir}/init_chimeras.lst') as f:
        init_list = f.read().splitlines()
    f.close
    
    ### Run multi-threads te-initiated function
    with Pool(threads) as pool:
        pool.starmap(init_chimeras, [(gene_id, itsct_mbed_TEs_up) for gene_id in init_list])

    if os.path.exists(f'{aln_dir}/TE-initiated_tmp.tsv'):
        if os.stat(f'{aln_dir}/TE-initiated_tmp.tsv').st_size > 0:
            os.rename(f'{aln_dir}/TE-initiated_tmp.tsv', f'{aln_dir}/TE-initiated-{group}.tsv')
            shutil.move(f'{aln_dir}/TE-initiated-{group}.tsv', f'{tmp}/TE-initiated-{group}.tsv')
            print(colored("Done!", "green", attrs=['bold']))
    else:
        print(colored('There are no TE-initiated transcripts!', "red"))
