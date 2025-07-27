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

def te_term(aln_dir, window):
    print("\n############################\n## TE-terminated analysis ##\n############################\n")

    clock = time()
    print(f'{clock}\tGenerating downstream region for expressed genes...')
    def create_down_bed(input_bed, output_bed, window):
        ### Function to create promoter region based on the windows size (3k default) | output = genes_UP_window.bed
        df = pd.read_csv(input_bed, sep="\t", header=None, names=['gene_id', 'start', 'end', 'id', 'dot', 'strand'])

        # Create new columns for promoter start and end
        df['down_end'] = df.apply(
            lambda row: max(1, row['end'] + window) if row['strand'] == '+' else row['start'],
            axis=1
        )
        df['down_start'] = df.apply(
            lambda row: row['end'] if row['strand'] == '+' else row['start'] - window,
            axis=1
        )

        # Reorder columns for BED format
        result = df[['gene_id', 'down_start', 'down_end', 'id', 'dot', 'strand']]

        # Save as BED
        result.to_csv(output_bed, sep="\t", header=False, index=False)
    ### Create bed file with promoter regions | output = genes_UP_window.bed
    create_down_bed(f'{aln_dir}/genes_total_expressed.bed', f'{aln_dir}/genes_DOWN_window.bed', window)
 
    ### Intersect promoters with TEs and get the list of the genes | output = init_chimeras.lst
    genes_TE_UP = intersection_any_bp_mbed(f'{aln_dir}/genes_DOWN_window.bed', f'{aln_dir}/chimeric_TEs.bed') 
    genes_TE_UP_IDs = get_IDs_from_bed(genes_TE_UP)
    genes_TE_UP_IDs["Name"].drop_duplicates().to_csv(f'{aln_dir}/term_chimeras.lst', header=None, index=False)
    
    print(colored("Done!", "green", attrs=['bold']))
    return genes_TE_UP

def term_chimeras(gene_id, itsct_mbed_TEs_up):
    from __main__ import aln_dir

    ### Extract all rows (TEs itersecting promoters)
    TEs_downstream = itsct_mbed_TEs_up.query('Name == @gene_id')[['Chromosome', 'Start_b', 'End_b', 'Name_b', 'Score_b', 'Strand_b']].drop_duplicates()
    TEs_downstream.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]

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
        for row in TEs_downstream.itertuples():
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
                cov = chim_reads_TE["Name"].isin(chim_reads_TE["Name"]).sum()
                if cov > 0:
                    chr_TE = row[1]
                    TE_strd = row[6]
                    TE_family = row[4]
                    strand_exon = str(exon_coord.iloc[0]["Strand"])

                    with open(f'{aln_dir}/TE-terminated_tmp.tsv', 'a') as te_init:
                        te_init.write(f'{gene_id}\t{strand_exon}\t{chr_gene}:{s_gene}-{e_gene}\t{TE_family}\t{TE_strd}\t{chr_TE}:{s_TE}-{e_TE}\t{cov}\n')
                    te_init.close()

def multicore_process_term(itsct_mbed_TEs_up, threads):
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import tmp

    clock = time()
    print(f'{clock}\tSearching for TE-terminated transcripts...')

    #TE-terminated chimeras
    with open(f'{aln_dir}/term_chimeras.lst') as f:
        init_list = f.read().splitlines()
    f.close
    
    ### Run multi-threads te-terminated function
    with Pool(threads) as pool:
        pool.starmap(term_chimeras, [(gene_id, itsct_mbed_TEs_up) for gene_id in init_list])

    if os.path.exists(f'{aln_dir}/TE-terminated_tmp.tsv'):
        if os.stat(f'{aln_dir}/TE-terminated_tmp.tsv').st_size > 0:
            os.rename(f'{aln_dir}/TE-terminated_tmp.tsv', f'{aln_dir}/TE-terminated-{group}.tsv')
            shutil.move(f'{aln_dir}/TE-terminated-{group}.tsv', f'{tmp}/TE-terminated-{group}.tsv')
            print(colored("Done!", "green", attrs=['bold']))
    else:
        print(colored('There are no TE-terminated transcripts!', "red"))



# from multiprocessing.dummy import Pool
# import concurrent.futures
# import argparse
# import datetime
# import subprocess
# from termcolor import colored
# import glob, os
# import pandas as pd
# import re
# import csv
# from io import StringIO
# from __main__ import *

# def te_term(aln_dir, group, out_group):
#     print("\n############################\n## TE-terminated analysis ##\n############################\n")
#     te_term=open(str(out_group + "/TE-terminated-" + str(group) + ".tsv"), 'w')

#     clock = time()
#     print(str(clock) + '\t' + "Generating downstream region for expressed genes...")
#     down_window = pd.read_table(str(aln_dir + '/genes_total_expressed.bed'), header=None, usecols=[0,1,2,3,4,5],names=['gene_id', 'start', 'end', 'id', 'dot', 'strand'])

#     with open(str(aln_dir + '/genes_DOWN_window.bed'), encoding='utf-8', mode='a') as down_window_file:
#         for index, row in down_window.iterrows():
#             scaf = repr(row['gene_id'])
#             start = repr(row['start'])
#             end = repr(row['end'])
#             part = repr(row['id'])
#             dot = repr(row['dot'])
#             strand = row['strand']
#             if strand == "+":
#                 num = (int(end) + int(args.window))
#                 down_window_file.write(scaf.replace("'","") + "\t" + str(end) + "\t" + str(num) + "\t" + part.replace("'","") + "\t" + dot.replace("'","") + "\t" + strand.replace("'","") + "\n")
#             else:
#                 num = (int(start) - int(args.window))
#                 if num < 1:
#                     num = 1
#                 down_window_file.write(scaf.replace("'", "") + "\t" + str(num) + "\t" + str(start) + "\t" + part.replace("'","") + "\t" + dot.replace("'", "") + "\t" + strand.replace("'", "") + "\n")
#     down_window_file.close

#     genes_DOWN_window = pybedtools.BedTool(str(aln_dir + '/genes_DOWN_window.bed'))
#     chimeric_TEs = pybedtools.BedTool(str(aln_dir + '/chimeric_TEs.bed'))
#     genes_TE_DOWN_mbed = chimeric_TEs.intersect(genes_DOWN_window, wa=True, wb=True, nonamecheck=True)
#     bed2string = str(genes_TE_DOWN_mbed)
#     bed2string = StringIO(bed2string)
#     global df
#     df = pd.read_csv(bed2string, sep="\t", header=None, usecols=[0,1,2,3,4,5,6,7,8,9,10,11],names=['scaf', 'start', 'end', 'ID', 'dot', 'strand','scaf2', 'start2', 'end2', 'ID2', 'dot2', 'strand2']).drop_duplicates()
#     df["ID2"].drop_duplicates().to_csv(str(aln_dir + '/term_chimeras.lst'), encoding='utf-8', header=None,index=False)
#     print(colored("Done!", "green", attrs=['bold']))

# def term_chimeras(term_list):
#     from __main__ import aln_dir
#     gene_reads_fwd = pybedtools.BedTool(aln_dir + '/gene_reads_fwd.bed')
#     TE_reads_fwd = pybedtools.BedTool(aln_dir + '/TE_reads_fwd.bed')
#     gene_reads_rev = pybedtools.BedTool(aln_dir + '/gene_reads_rev.bed')
#     TE_reads_rev = pybedtools.BedTool(aln_dir + '/TE_reads_rev.bed')

#     exon = pd.read_csv(str(aln_dir + '/chim_exons.bed'), sep="\t", usecols=[0,1,2,3,4,5], names=['scaf_exon', 'start_exon', 'end_exon', 'ID_exon', 'dot_exon', 'strand_exon'])

#     gene_id = term_list
#     TEs_upstream = df.query('ID2 == @gene_id')
#     gene_coord = df.query('ID2 == @gene_id').iloc[:, 6:12].drop_duplicates()
#     chr_gene = gene_coord['scaf2'].to_string(index=False).replace(" ","")
#     s_gene = gene_coord['start2'].to_string(index=False).replace(" ","")
#     e_gene = gene_coord['end2'].to_string(index=False).replace(" ","")
#     strand = gene_coord['strand2'].to_string(index=False).replace(" ","")
#     exon_coord = exon.query('ID_exon == @gene_id').to_csv(sep='\t', encoding='utf-8', header=None,index=False)
#     exon_bed = pybedtools.BedTool(str(exon_coord), from_string=True)

#     if exon_coord:
#         exons_coord = StringIO(exon_coord)
#         strand_exon = pd.read_csv(exons_coord, sep = "\t", header=None, usecols=[5]).iloc[0].item()

#         for row in TEs_upstream.itertuples():
#             TE_tmp = row[1:7]
#             TE_tmp = pd.DataFrame(TE_tmp).transpose()

#             chr_TE = row[1]
#             s_TE = row[2]
#             e_TE = row[3]
#             TE_family = row[4]
#             dot = row[5]
#             TE_strd = row[6]

#             merging_TE_reads = None
#             merging_gene_reads = None
#             if strand == "+":
#                 if int(s_TE) < int(s_gene):
#                     s_TE = int(s_gene)
#                 TE_coord = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
#                 TE_bed = pybedtools.BedTool(TE_coord, from_string=True)
#                 read_genes = gene_reads_fwd.intersect(exon_bed, wa=True, wb=True, nonamecheck=True)
#                 intersect_gene = str(read_genes)
#                 intersect_gene = StringIO(intersect_gene)
#                 reads_gene_col = pd.read_table(intersect_gene, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

#                 reads_TE = TE_reads_fwd.intersect(TE_bed, wa=True, wb=True, nonamecheck=True, f=str(args.overlap))
#                 intersect_TE = str(reads_TE)
#                 intersect_TE = StringIO(intersect_TE)
#                 reads_TE_col = pd.read_table(intersect_TE, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()
#             else:
#                 if int(e_TE) > int(e_gene):
#                     e_TE = int(e_gene)
#                 TE_coord = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
#                 TE_bed = pybedtools.BedTool(TE_coord, from_string=True)
#                 read_genes = gene_reads_rev.intersect(exon_bed, wa=True, wb=True, nonamecheck=True)
#                 intersect_gene = str(read_genes)
#                 intersect_gene = StringIO(intersect_gene)
#                 reads_gene_col = pd.read_table(intersect_gene, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

#                 reads_TE = TE_reads_rev.intersect(TE_bed, wa=True, wb=True, nonamecheck=True, f=str(args.overlap))
#                 intersect_TE = str(reads_TE)
#                 intersect_TE = StringIO(intersect_TE)
#                 reads_TE_col = pd.read_table(intersect_TE, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

#             if not reads_gene_col.empty:
#                 if merging_gene_reads is not None:
#                     merging_gene_reads = merging_gene_reads.append(reads_gene_col, ignore_index=True).drop_duplicates()
#                 else:
#                     merging_gene_reads = reads_gene_col

#             if not reads_TE_col.empty:
#                 if merging_TE_reads is not None:
#                     merging_TE_reads = merging_TE_reads.append(reads_TE_col, ignore_index=True).drop_duplicates()
#                 else:
#                     merging_TE_reads = reads_TE_col

#             if merging_TE_reads is not None:
#                 if merging_gene_reads is not None:
#                     list_TE_reads = merging_TE_reads['read_ID'].tolist()
#                     match = merging_gene_reads[merging_gene_reads['read_ID'].isin(list_TE_reads)].drop_duplicates()
#                     cov = len(match)
#                     if cov > 0:
#                         with open(str("TE-terminated_tmp.tsv"), "a") as te_term:
#                             print(str(gene_id) + "\t" + str(strand_exon) + "\t" + str(chr_gene) + ":" + str(s_gene) + "-" + str(e_gene) + "\t" + str(TE_family) + "\t" + str(TE_strd) + "\t" + str(chr_TE) + ":" + str(s_TE) + "-" + str(e_TE) + "\t" + str(cov), file=te_term)
#                         te_term.close()

# def multicore_process_term():
#     from __main__ import aln_dir
#     from __main__ import group
#     from __main__ import out_group

#     clock = time()
#     print(str(clock) + '\t' + "Searching for TE-terminated transcripts...")

#     #TE-terminated chimeras
#     with open(str(aln_dir + '/term_chimeras.lst')) as f:
#         term_list = f.read().splitlines()
#     f.close
#     pool = Pool(processes=int(args.threads))
#     pool.map(term_chimeras, term_list)
#     pool.close()

#     term_output = str('TE-terminated_tmp.tsv')
#     if os.path.exists(term_output) == True:
#         if os.stat(str('TE-terminated_tmp.tsv')).st_size > 0:
#             os.rename('TE-terminated_tmp.tsv', str('TE-terminated-' + str(group) + '.tsv'))
#             copy(str('TE-terminated-' + str(group) + '.tsv'), out_group)
#             os.remove('TE-terminated-' + str(group) + '.tsv')
#             print(colored("Done!", "green", attrs=['bold']))
#             pybedtools.cleanup(remove_all=True)
#     else:
#         print(colored('There are no TE-terminated transcripts!', "red"))
