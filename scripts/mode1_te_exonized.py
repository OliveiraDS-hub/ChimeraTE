import multiprocessing
import concurrent.futures
from multiprocessing.dummy import Pool
import argparse
import datetime
from termcolor import colored
import subprocess
import glob, os
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *
import pybedtools
from multiprocessing import Pool
import os.path

def te_exon_embedded(aln_dir, group, out_group):
    print("\n###########################\n## TE-exonized analysis ##\n###########################\n")

    ##Identify TEs within genes
    clock = time()
    print(str(clock) + '\t' + "Identifying embedded, overlapped, and intronic TEs...")

    intersection(str(aln_dir + '/chimeric_TEs.bed'), str(aln_dir + '/genes_total_expressed.bed'), str(aln_dir + '/genes_TE_INSIDE.mbed'), str('0.1')) #All genes with TEs >= 10% inside
    genes_TE_INSIDE = pd.read_csv(str(aln_dir + '/genes_TE_INSIDE.mbed'), sep="\t").iloc[:,9].drop_duplicates().tolist() #Gene IDs
    overlap(str(aln_dir + '/chim_exons.bed'), genes_TE_INSIDE, str(aln_dir + '/gene_TE_INSIDE-ALL_exons.bed')) #All exons from genes with TEs >= 10% inside

    if os.stat(str(aln_dir + '/gene_TE_INSIDE-ALL_exons.bed')).st_size > 0:
        intersection(str(aln_dir + '/chimeric_TEs.bed'), str(aln_dir + '/gene_TE_INSIDE-ALL_exons.bed'), str(aln_dir + '/gene_TE_INSIDE-exons.mbed'), str('0.1'))

        ##Embedded
        pd.read_csv(str(aln_dir + '/genes_TE_INSIDE.mbed'), sep = "\t").iloc[:,0:6].drop_duplicates().to_csv(aln_dir + '/TEs_inside_genes.bed', sep='\t', header=None, index=False) #All TEs with >= 10% inside genes
        pd.read_csv(str(aln_dir + '/gene_TE_INSIDE-exons.mbed'), sep = "\t").iloc[:,6:12].drop_duplicates().to_csv(aln_dir + '/gene_TE_INSIDE-exons.bed', sep='\t', header=None, index=False) #All TEs with >= 10% inside genes

        intersection(str(aln_dir + '/TEs_inside_genes.bed'), str(aln_dir + '/gene_TE_INSIDE-exons.bed'), str(aln_dir + '/embedded_exons.bed'), str('1')) # TEs and exons with TEs = 100% embedded
        global genes_TE_embedded
        genes_TE_embedded = pd.read_csv(str(aln_dir + '/embedded_exons.bed'), sep="\t").iloc[:,9].drop_duplicates().tolist() #List of genes with embedded TEs on exons

        global genes_TE_INSIDE_bed
        genes_TE_INSIDE_bed = pd.read_csv(str(aln_dir + '/genes_TE_INSIDE.mbed'), sep = "\t", usecols=[0,1,2,3,4,5,6,7,8,9,10,11],names=['scaf_TE', 'start_TE', 'end_TE', 'ID_TE', 'dot_TE', 'strand_TE','scaf_gene', 'start_gene', 'end_gene', 'ID_gene', 'dot_gene', 'strand_gene'])
        genes_TE_INSIDE_bed = genes_TE_INSIDE_bed[genes_TE_INSIDE_bed['ID_gene'].isin(genes_TE_embedded)].drop_duplicates()
        genes_TE_INSIDE_bed["ID_gene"].drop_duplicates().to_csv(str(aln_dir + '/exon_emb.lst'), encoding='utf-8', header=None,index=False)

        global embedded_exons
        embedded_exons = pd.read_csv(str(aln_dir + '/embedded_exons.bed'), sep = "\t", usecols=[6,7,8,9,10,11],names=['scaf_exon', 'start_exon', 'end_exon', 'ID_exon', 'dot_exon', 'strand_exon']).drop_duplicates()
        print(colored("Done!", "green", attrs=['bold']))    # OK    # OK (old-school method) # OK (old-school method)

def embedded_mp(all_emb_IDs):
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import out_group
    gene_reads_fwd = pybedtools.BedTool(aln_dir + '/gene_reads_fwd.bed')
    TE_reads_fwd = pybedtools.BedTool(aln_dir + '/TE_reads_fwd.bed')
    gene_reads_rev = pybedtools.BedTool(aln_dir + '/gene_reads_rev.bed')
    TE_reads_rev = pybedtools.BedTool(aln_dir + '/TE_reads_rev.bed')

    gene_id = all_emb_IDs

    gene_coord = genes_TE_INSIDE_bed.query('ID_gene == @gene_id').iloc[:, 6:12].drop_duplicates()
    chr_gene = gene_coord['scaf_gene'].to_string(index=False).replace(" ","")
    s_gene = gene_coord['start_gene'].to_string(index=False).replace(" ","")
    e_gene = gene_coord['end_gene'].to_string(index=False).replace(" ","")
    strand_gene = gene_coord['strand_gene'].to_string(index=False).replace(" ","")

    exon_coord = embedded_exons.query('ID_exon == @gene_id').to_csv(sep='\t', encoding='utf-8', header=None,index=False)
    exon_bed = pybedtools.BedTool(str(exon_coord), from_string=True)

    TE_embedded = genes_TE_INSIDE_bed.query('ID_gene == @gene_id').iloc[:, 0:6].drop_duplicates()

    for row in TE_embedded.itertuples():
        chr_TE = row[1]
        s_TE = row[2]
        e_TE = row[3]
        TE_family = row[4]
        dot = row[5]
        TE_strd = row[6]

        TE_coord = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
        TE_bed = pybedtools.BedTool(str(TE_coord), from_string=True)

        merging_TE_reads = None
        merging_gene_reads = None

        if strand_gene == "+":
            read_genes = gene_reads_fwd.intersect(exon_bed, wa=True, wb=True, nonamecheck=True)
            intersect_gene = str(read_genes)
            intersect_gene = StringIO(intersect_gene)
            reads_gene_col = pd.read_table(intersect_gene, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

            reads_TE = TE_reads_fwd.intersect(TE_bed, wa=True, wb=True, nonamecheck=True, f=str(args.overlap))
            intersect_TE = str(reads_TE)
            intersect_TE = StringIO(intersect_TE)
            reads_TE_col = pd.read_table(intersect_TE, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

        else:
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
                    with open(str(out_group + '/TE-exonized-' + str(group) + '.tsv'), "a") as te_exon:
                        print(str(gene_id) + "\t" + str(strand_gene) + "\t" + str(chr_gene) + ":" + str(s_gene) + "-" + str(e_gene) + "\t" + str(TE_family) + "\t" + str(TE_strd) + "\t" + str(chr_TE) + ":" + str(s_TE) + "-" + str(e_TE) + "\t" + str(cov) + "\t" + str("Embedded"), file=te_exon)
                    te_exon.close()


def prep_overlapped(aln_dir, group, out_group):

    ##Identify TEs overlapping exons
    exons_wsome_overlapping = pd.read_csv(str(aln_dir + '/gene_TE_INSIDE-exons.bed'), sep ='\t', usecols=[0,1,2,3,4,5],names=['scaf_exon', 'start_exon', 'end_exon', 'ID_exon', 'dot_exon', 'strand_exon']).drop_duplicates() #All exons from genes with TEs >= 10% inside

    pd.concat([exons_wsome_overlapping,embedded_exons]).drop_duplicates(keep=False).to_csv(aln_dir + '/non_embedded_exons.bed', sep = "\t", header=None, index=False)

    TEs_embedded = pd.read_csv(str(aln_dir + '/genes_TE_INSIDE.mbed'), sep ='\t', usecols=[0,1,2,3,4,5],names=['scaf_TE', 'start_TE', 'end_TE', 'ID_TE', 'dot_TE', 'strand_TE']).drop_duplicates()
    TEs_inside = pd.read_csv(str(aln_dir + '/embedded_exons.bed'), sep = "\t", usecols=[0,1,2,3,4,5],names=['scaf_TE', 'start_TE', 'end_TE', 'ID_TE', 'dot_TE', 'strand_TE']).drop_duplicates()#.to_csv(aln_dir + '/gene_TE_INSIDE-exons.bed', sep='\t', header=None, index=False)
    pd.concat([TEs_inside,TEs_embedded]).drop_duplicates(keep=False).to_csv(aln_dir + '/TEs_overlapping.bed', sep = "\t", header=None, index=False)


    intersection(str(aln_dir + '/non_embedded_exons.bed'), str(aln_dir + '/TEs_overlapping.bed'), str(aln_dir + '/gene_TE_OVERLAPPED-exons.mbed')) #All TEs and exons from genes with TEs overlapping
    global genes_TE_overlapping
    if os.stat(str(aln_dir + '/gene_TE_OVERLAPPED-exons.mbed')).st_size > 0:
        genes_TE_overlapping = pd.read_csv(str(aln_dir + '/gene_TE_OVERLAPPED-exons.mbed'), sep = '\t', usecols=[0,1,2,3,4,5,6,7,8,9,10,11],names=['scaf_exon', 'start_exon', 'end_exon', 'ID_exon', 'dot_exon', 'strand_exon', 'scaf_TE', 'start_TE', 'end_TE', 'ID_TE', 'dot_TE', 'strand_TE'])
        global exons_overlapping
        exons_overlapping = pd.read_csv(str(aln_dir + '/gene_TE_OVERLAPPED-exons.mbed'), sep = '\t', usecols=[0,1,2,3,4,5],names=['scaf_exon', 'start_exon', 'end_exon', 'ID_exon', 'dot_exon', 'strand_exon'])

        pd.read_csv(str(aln_dir + '/gene_TE_OVERLAPPED-exons.mbed'), sep = "\t").iloc[:,6:12].drop_duplicates().to_csv(aln_dir + '/TEs_overlapped_exons.bed', sep='\t', header=None, index=False) #All TEs with >= 10% inside genes
        intersection(str(aln_dir + '/TEs_overlapped_exons.bed'), str(aln_dir + '/genes_total_expressed.bed'), str(aln_dir + '/genes_TE_overlapping.mbed'), str('0.01')) #All TEs and exons from genes with TEs overlapping
        genes_TE_overlapping = pd.read_csv(str(aln_dir + '/genes_TE_overlapping.mbed'), sep = '\t', usecols=[0,1,2,3,4,5,6,7,8,9,10,11],names=['scaf_TE', 'start_TE', 'end_TE', 'ID_TE', 'dot_TE', 'strand_TE', 'scaf_gene', 'start_gene', 'end_gene', 'ID_gene', 'dot_gene', 'strand_gene'])
        genes_TE_overlapping["ID_gene"].drop_duplicates().to_csv(str(aln_dir + '/exon_overlapped.lst'), encoding='utf-8', header=None,index=False)

def overlapped_mp(all_overlap_IDs):
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import out_group
    gene_reads_fwd = pybedtools.BedTool(aln_dir + '/gene_reads_fwd.bed')
    TE_reads_fwd = pybedtools.BedTool(aln_dir + '/TE_reads_fwd.bed')
    gene_reads_rev = pybedtools.BedTool(aln_dir + '/gene_reads_rev.bed')
    TE_reads_rev = pybedtools.BedTool(aln_dir + '/TE_reads_rev.bed')

    gene_id = all_overlap_IDs

    gene_coord = genes_TE_overlapping.query('ID_gene == @gene_id').iloc[:, 6:12].drop_duplicates()
    chr_gene = gene_coord['scaf_gene'].to_string(index=False).replace(" ","")
    s_gene = gene_coord['start_gene'].to_string(index=False).replace(" ","")
    e_gene = gene_coord['end_gene'].to_string(index=False).replace(" ","")
    strand_gene = gene_coord['strand_gene'].to_string(index=False).replace(" ","")

    exon_coord = exons_overlapping.query('ID_exon == @gene_id').to_csv(sep='\t', encoding='utf-8', header=None,index=False)
    exon_bed = pybedtools.BedTool(str(exon_coord), from_string=True)

    TEs_overlapping = genes_TE_overlapping.query('ID_gene == @gene_id').iloc[:, 0:6].drop_duplicates()

    for row in TEs_overlapping.itertuples():
        chr_TE = row[1]
        s_TE = row[2]
        e_TE = row[3]
        TE_family = row[4]
        dot = row[5]
        TE_strd = row[6]

        TE_coord = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
        TE_bed = pybedtools.BedTool(str(TE_coord), from_string=True)

        merging_TE_reads = None
        merging_gene_reads = None
        if int(s_TE) < int(s_gene):
            s_TE = s_gene
        if int(e_TE) > int(e_gene):
            e_TE = e_gene

        exon_coord = exons_overlapping.query('ID_exon == @gene_id')

        for row in exon_coord.itertuples():
            chr_exon = row[1]
            s_exon = row[2]
            e_exon = row[3]

            if int(s_TE) >= int(s_exon) and int(e_TE) <= int(e_exon):
                continue
            else:
                if int(s_TE) <= int(e_exon) and int(s_TE) >= int(s_exon) and int(e_TE) > int(e_exon):
                    s_TE = e_exon
                    TE_coord = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
                    TE_bed = pybedtools.BedTool(str(TE_coord), from_string=True)
                elif int(s_TE) < int(s_exon) and int(e_TE) >= int(s_exon) and int(e_TE) <= int(e_exon):
                    e_TE = s_exon
                    TE_coord = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
                    TE_bed = pybedtools.BedTool(str(TE_coord), from_string=True)
                elif int(s_TE) < int(s_exon) and int(e_TE) >= int(e_exon):
                    e1_TE = s_exon; s2_TE = e_exon
                    TE_coord1 = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e1_TE) + "\t" + str(TE_family)
                    TE_coord2 = str(chr_TE) + "\t" + str(s2_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
                    TE_coord_merged = str(TE_coord1) + '\n' + str(TE_coord2)
                    TE_bed = pybedtools.BedTool(str(TE_coord_merged), from_string=True)

            if TE_bed is not None:
                if strand_gene == "+":
                    read_genes = gene_reads_fwd.intersect(exon_bed, wa=True, wb=True, nonamecheck=True)
                    intersect_gene = str(read_genes)
                    intersect_gene = StringIO(intersect_gene)
                    reads_gene_col = pd.read_table(intersect_gene, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

                    reads_TE = TE_reads_fwd.intersect(TE_bed, wa=True, wb=True, nonamecheck=True, f=str(args.overlap))
                    intersect_TE = str(reads_TE)
                    intersect_TE = StringIO(intersect_TE)
                    reads_TE_col = pd.read_table(intersect_TE, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()
                else:
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
                    with open(str(out_group + '/TE-exonized-' + str(group) + '.tsv'), "a") as te_exon:
                        print(str(gene_id) + "\t" + str(strand_gene) + "\t" + str(chr_gene) + ":" + str(s_gene) + "-" + str(e_gene) + "\t" + str(TE_family) + "\t" + str(TE_strd) + "\t" + str(chr_TE) + ":" + str(s_TE) + "-" + str(e_TE) + "\t" + str(cov) + "\t" + str("Overlapped"), file=te_exon)
                    te_exon.close()

def prep_intronic(aln_dir, group, out_group):
    from __main__ import aln_dir

    ##Identify TEs in intronic region
    intersection(str(aln_dir + '/chimeric_TEs.bed'), str(aln_dir + '/genes_total_expressed.bed'), str(aln_dir + '/genes_TE_INSIDE_any_overlap.mbed')) #All genes with TEs >= 10% inside

    pd.read_csv(str(aln_dir + '/genes_TE_INSIDE_any_overlap.mbed'), sep = "\t").iloc[:,0:6].drop_duplicates().to_csv(aln_dir + '/TEs_inside_genes_any_overlap.bed', sep='\t', header=None, index=False) #All TEs with >= 10% inside genes

    all_TEs_within_genes = pd.read_csv(str(aln_dir + '/TEs_inside_genes_any_overlap.bed'), sep = "\t", usecols=[0,1,2,3,4,5],names=['scaf_TE', 'start_TE', 'end_TE', 'ID_TE', 'dot_TE', 'strand_TE']).drop_duplicates()
    pd.read_csv(str(aln_dir + '/genes_TE_INSIDE.mbed'), sep = "\t").iloc[:,6:12].drop_duplicates().to_csv(aln_dir + '/genes_with_TEs.bed', sep='\t', header=None, index=False)

    intersection(str(aln_dir + '/chimeric_TEs.bed'), str(aln_dir + '/gene_TE_INSIDE-ALL_exons.bed'), str(aln_dir + '/gene_TE_overlapping-exons.mbed'))
    overlapped_TEs = pd.read_csv(str(aln_dir + '/gene_TE_overlapping-exons.mbed'), sep = "\t", usecols=[0,1,2,3,4,5],names=['scaf_TE', 'start_TE', 'end_TE', 'ID_TE', 'dot_TE', 'strand_TE']).drop_duplicates()

    pd.concat([all_TEs_within_genes,overlapped_TEs]).drop_duplicates(keep=False).to_csv(aln_dir + '/intronic_TEs.bed', sep = "\t", header=None, index=False)

    intersection(str(aln_dir + '/genes_with_TEs.bed'), str(aln_dir + '/intronic_TEs.bed'), str(aln_dir + '/gene_intronic_boundaries_TEs.mbed'))
    global TEs_intron_nearby
    TEs_intron_nearby = pd.read_csv(str(aln_dir + '/gene_intronic_boundaries_TEs.mbed'), sep = "\t", usecols=[0,1,2,3,4,5,6,7,8,9,10,11],names=['scaf_gene', 'start_gene', 'end_gene', 'ID_gene', 'dot_gene', 'strand_gene', 'scaf_TE', 'start_TE', 'end_TE', 'ID_TE', 'dot_TE', 'strand_TE'])

    genes_intronic_TEs = pd.read_csv(str(aln_dir + '/gene_intronic_boundaries_TEs.mbed'), sep="\t").iloc[:,3].drop_duplicates().tolist() #Gene IDs
    overlap(str(aln_dir + '/chim_exons.bed'), genes_intronic_TEs, str(aln_dir + '/gene_TE_intronic_exons.bed'))

    global gene_TE_intronic_exons
    gene_TE_intronic_exons = pd.read_csv(str(aln_dir + '/gene_TE_intronic_exons.bed'), sep = "\t", usecols=[0,1,2,3,4,5],names=['scaf_exon', 'start_exon', 'end_exon', 'ID_exon', 'dot_exon', 'strand_exon']).drop_duplicates()

    gene_TE_intronic_exons["ID_exon"].drop_duplicates().to_csv(aln_dir + '/exon_intron_list.lst', header=None, index=False) # OK

def intronic_mp(all_gene_IDs):
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import out_group
    gene_reads_fwd = pybedtools.BedTool(aln_dir + '/gene_reads_fwd.bed')
    TE_reads_fwd = pybedtools.BedTool(aln_dir + '/TE_reads_fwd.bed')
    gene_reads_rev = pybedtools.BedTool(aln_dir + '/gene_reads_rev.bed')
    TE_reads_rev = pybedtools.BedTool(aln_dir + '/TE_reads_rev.bed')

    gene = all_gene_IDs
    TEs_intron_nearby_specific_gene = TEs_intron_nearby.query('ID_gene == @gene')
    for row in TEs_intron_nearby_specific_gene.itertuples():
        chr_gene = row[1]
        s_gene = row[2]
        e_gene = row[3]
        gene_ID = row[4]
        dot_gene = row[5]
        gene_strd = row[6]
        chr_TE = row[7]
        s_TE = row[8]
        e_TE = row[9]
        dot_TE = row[10]
        TE_family = row[10]
        TE_strd = row[12]

        TE_coord = str(chr_TE) + "\t" + str(s_TE) + "\t" + str(e_TE) + "\t" + str(TE_family)
        TE_bed = pybedtools.BedTool(TE_coord, from_string=True)

        exons_coord = gene_TE_intronic_exons.query('ID_exon == @gene_ID').to_csv(index=False, header=False, sep='\t')
        exons_bed = pybedtools.BedTool(exons_coord, from_string=True)

        if exons_bed:
            exons_coord = StringIO(exons_coord)
            strand_exon = pd.read_csv(exons_coord, sep = "\t", header=None, usecols=[5]).iloc[0].item()
        else:
            continue

        if TE_bed:
            merging_TE_reads = None
            merging_gene_reads = None

            if strand_exon == "+":
                read_genes = gene_reads_fwd.intersect(exons_bed, wa=True, wb=True, nonamecheck=True)
                intersect_gene = str(read_genes)
                intersect_gene = StringIO(intersect_gene)
                reads_gene_col = pd.read_table(intersect_gene, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()

                reads_TE = TE_reads_fwd.intersect(TE_bed, wa=True, wb=True, nonamecheck=True, f=str(args.overlap))
                intersect_TE = str(reads_TE)
                intersect_TE = StringIO(intersect_TE)
                reads_TE_col = pd.read_table(intersect_TE, sep="\t", header=None, usecols=[3],names=['read_ID']).replace('/1', '', regex=True).replace('/2', '', regex=True).drop_duplicates()
            else:
                read_genes = gene_reads_rev.intersect(exons_bed, wa=True, wb=True, nonamecheck=True)
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

        if merging_TE_reads is None:
                continue
        else:
            if merging_gene_reads is not None:
                list_TE_reads = merging_TE_reads['read_ID'].tolist()
                match = merging_gene_reads[merging_gene_reads['read_ID'].isin(list_TE_reads)].drop_duplicates()
                cov = len(match)
                if cov > 0:
                    with open(str(out_group + '/TE-exonized-' + str(group) + '.tsv'), "a") as te_exon:
                        print(str(gene_ID) + "\t" + str(strand_exon) + "\t" + str(chr_gene) + ":" + str(s_gene) + "-" + str(e_gene) + "\t" + str(TE_family) + "\t" + str(TE_strd) + "\t" + str(chr_TE) + ":" + str(s_TE) + "-" + str(e_TE) + "\t" + str(cov) + "\t" + str("Intronic"), file=te_exon)
                    te_exon.close()

def multicore_process_exon():
    from __main__ import aln_dir
    from __main__ import group
    from __main__ import out_group

    ### TE-exonized EMBEDDED chimeras
    clock = time()
    print(str(clock) + "\t" + "Searching Embedded-TE chimeras...")
    list = str(aln_dir + '/exon_emb.lst')

    if os.path.exists(list) == True:
        with open(str(aln_dir + '/exon_emb.lst')) as f:
            all_emb_IDs = f.read().splitlines()
        f.close
        pool = Pool(processes=int(args.threads))
        pool.map(embedded_mp, all_emb_IDs)
        pool.close()
        print(colored("Done!", "green", attrs=['bold']))
        pybedtools.cleanup(remove_all=True)
    else:
        print(colored('There are no TEs embedded into exons with chimeric reads!', "red"))


    ### TE-exonized OVERLAPPED chimeras
    clock = time()
    print(str(clock) + '\t' + "Searching Overlapped-TE chimeras...")
    list = str(aln_dir + '/exon_overlapped.lst')

    if os.path.exists(list) == True:
        with open(str(aln_dir + '/exon_overlapped.lst')) as f:
            all_overlap_IDs = f.read().splitlines()
        f.close
        pool = Pool(processes=int(args.threads))
        pool.map(overlapped_mp, all_overlap_IDs)
        pool.close()
        print(colored("Done!", "green", attrs=['bold']))
        pybedtools.cleanup(remove_all=True)
    else:
        print(colored('There are no TEs overlapping exons with chimeric reads!', "red"))

    ### TE-exonized INTRONIC chimeras
    clock = time()
    print(str(clock) + '\t' + "Searching Intronic-TE chimeras...")
    list = str(aln_dir + '/exon_intron_list.lst')

    if os.path.exists(list) == True:
        with open(str(aln_dir + '/exon_intron_list.lst')) as f:
            all_gene_IDs = f.read().splitlines()
        f.close
        pool = Pool(processes=int(args.threads))
        pool.map(intronic_mp, all_gene_IDs)
        pool.close()
        print(colored("Done!", "green", attrs=['bold']))
        pybedtools.cleanup(remove_all=True)
    else:
        print(colored('There are no TEs into introns with chimeric reads!', "red"))
