import argparse
import datetime
import subprocess
import glob
import pandas as pd
import re
import csv
from io import StringIO
from __main__ import *
import pybedtools
import contextlib
import sys
import os


def alignment_func(out_dir,group,aln_dir,mate1,mate2):
    print(f"Running analysis with ------------------------------------------> {group}\n")
    clock = time()
    print(f"{clock}\tPerforming alignment")
    star_threads = int(int(args.threads) * 0.8)
    if star_threads % 2 != 0:
        star_threads -= 1

    ### Perform STAR alignment
    if os.path.exists(f"{aln_dir}/{group}_Aligned.sortedByCoord.out.bam"):
        print(f"{group} bam file from STAR has been found at: {aln_dir}/{group}_Aligned.sortedByCoord.out.bam\n"
              f"skipping alignment...")
    elif str(args.index) != "None":
        subprocess.call(['STAR', '--genomeDir', str(args.index), '--runThreadN', str(star_threads), "--readFilesCommand", str("zcat"), \
        '--readFilesIn', str(mate1), str(mate2), '--outSAMtype', str("BAM"), str("SortedByCoordinate"), "--outFileNamePrefix", str(f"{aln_dir}/{group}_")], stdout=subprocess.DEVNULL)
    else:
        subprocess.call(['STAR', '--genomeDir', str(f"{out_dir}/index"), '--runThreadN', str(star_threads), "--readFilesCommand", str("zcat"), \
        '--readFilesIn', str(mate1), str(mate2), '--outSAMtype', str("BAM"), str("SortedByCoordinate"), "--outFileNamePrefix", str(f"{aln_dir}/{group}_")], stdout=subprocess.DEVNULL)

    ### Get unique aligned reads
    subprocess.call(['samtools', 'view', '-@', str(args.threads), '-b', '-q', '255', str(f"{aln_dir}/{group}_Aligned.sortedByCoord.out.bam")], stdout=open(str(f"{aln_dir}/accepted_hits.bam"), 'w'))

    ### Convert bam to bed
    subprocess.call(['bedtools', 'bamtobed', '-split', '-i', str(f"{aln_dir}/accepted_hits.bam")], stdout=open(str(f"{aln_dir}/accepted_hits.bed"), 'w'))

    ### Get reads aligned on the forward strand
    subprocess.call(['samtools', 'view', '-@', str(args.threads), '-b', '-f', '128', '-F', '16', str(f"{aln_dir}/accepted_hits.bam")], stdout=open(str(f"{aln_dir}/fwd1_f.bam"), 'w'))
    subprocess.call(['samtools', 'view', '-@', str(args.threads), '-b', '-f', '80', str(f"{aln_dir}/accepted_hits.bam")], stdout=open(str(f"{aln_dir}/fwd2_f.bam"), 'w'))

    ### Combine alignments that originate on the forward strand.
    subprocess.call(['samtools', 'merge', '-@', str(args.threads), '-f', str(f"{aln_dir}/fwd.bam"), str(f"{aln_dir}/fwd1_f.bam"), str(f"{aln_dir}/fwd2_f.bam")])


    ### Reverse strand
    subprocess.call(['samtools', 'view', '-@', str(args.threads), '-b', '-f', '144', str(f"{aln_dir}/accepted_hits.bam")], stdout=open(str(f"{aln_dir}/rev1_r.bam"), 'w'))
    subprocess.call(['samtools', 'view', '-@', str(args.threads), '-b', '-f', '64', '-F', '16', str(f"{aln_dir}/accepted_hits.bam")], stdout=open(str(f"{aln_dir}/rev2_r.bam"), 'w'))

    ### Combine alignments that originate on the reverse strand.
    subprocess.call(['samtools', 'merge', '-@', str(args.threads), '-f', str(f"{aln_dir}/rev.bam"), str(f"{aln_dir}/rev1_r.bam"), str(f"{aln_dir}/rev2_r.bam")])

    ### Foward bam to bed
    subprocess.call(['bedtools', 'bamtobed', '-split', '-i', str(f"{aln_dir}/fwd.bam")], stdout=open(str(f"{aln_dir}/fwd.bed"), 'w'))
    fwd_bed = pybedtools.BedTool(str(f"{aln_dir}/fwd.bed"))

    ### Reverse bam to bed
    subprocess.call(['bedtools', 'bamtobed', '-split', '-i', str(f"{aln_dir}/rev.bam")], stdout= open(str(f"{aln_dir}/rev.bed"), 'w'))
    rev_bed = pybedtools.BedTool(str(aln_dir + '/rev.bed'))
    print(colored("Done!", "green", attrs=['bold']))


    clock = time()
    print(f"{clock}\tGenes expression")
    fpkm_dir = create_dir(str(out_dir + '/' + group + '/alignment/fpkm'))

    if str(args.strand) == str("rf-stranded"):
        subprocess.call(['cufflinks', str(f"{aln_dir}/accepted_hits.bam"), '-p', str(args.threads), '-G', str(f"{tmp}/gtf_file.gtf"), '-o', str(fpkm_dir), '--quiet', '--library-type', 'fr-firststrand'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.call(['cufflinks', str(f"{aln_dir}/accepted_hits.bam"), '-p', str(args.threads), '-G', str(f"{tmp}/gtf_file.gtf"), '-o', str(fpkm_dir), '--quiet', '--library-type', 'fr-secondstrand'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    TEfile = pybedtools.BedTool(str(tmp + '/TE_file.bed'))
    acc_hits = pybedtools.BedTool(str(aln_dir + '/accepted_hits.bed'))
    exp_TEs = TEfile.intersect(acc_hits, wa=True, nonamecheck=True)
    exp_TEs = dropdup_bed(exp_TEs)
    with open(str(aln_dir + '/expressed_TEs.bed'), "w") as f:
        print(str(exp_TEs), file = f)

    expressed_TEs = pybedtools_intersection(str(tmp + '/TE_file.bed'), str(aln_dir + '/accepted_hits.bed'))
    expressed_TEs = dropdup_bed(expressed_TEs)
    expressed_TEs = pybedtools.BedTool(expressed_TEs, from_string=True)
    print(colored("Done!", "green", attrs=['bold']))

    clock = time()
    print(f"{clock}\tGetting fpkm")
    fpkm_gene = pd.read_csv(str(f"{out_dir}/{group}/alignment/fpkm/genes.fpkm_tracking"), sep="\t", usecols=[0,11])
    fpkm_gene['tracking_id'] = fpkm_gene['tracking_id'].str.replace("gene-", '')
    fpkm_gene.loc[fpkm_gene['FPKM_conf_hi'] > int(args.fpkm)].loc[:,['tracking_id']].drop_duplicates().to_csv(f"{aln_dir}/genes_expressed_IDs.lst", header=None, index=False)

    with open(str(aln_dir + '/genes_expressed_IDs.lst')) as f:
        list_genes = f.read().splitlines()

    overlap(str(tmp + '/gene_coord.bed'), list_genes, str(f"{aln_dir}/genes_total_expressed.bed"))
    intersected_genes = pybedtools_intersection(str(f"{aln_dir}/accepted_hits.bed"), str(f"{tmp}/gene_coord.bed"))

    print(colored("Done!", "green", attrs=['bold']))

    clock = time()
    print(f"{clock}\tStrand-specific expression analysis")

    ### TE reads
    TE_reads_fwd_all = fwd_bed.intersect(expressed_TEs, wa=True, nonamecheck=True)
    TE_reads_rev_all = rev_bed.intersect(expressed_TEs, wa=True, nonamecheck=True)

    fwd_TE = get_IDs_from_bed(TE_reads_fwd_all)
    rev_TE = get_IDs_from_bed(TE_reads_rev_all)

    merged_TE = pd.concat([fwd_TE, rev_TE])
    merged_TE["ID"].str.replace('/1', '').str.replace('/2', '').drop_duplicates().to_csv(f"{aln_dir}/TE_reads.lst", header=None, index=False)

    ### Gene reads
    genes_total_expressed = pybedtools.BedTool(str(f"{aln_dir}/genes_total_expressed.bed"))
    gene_reads_fwd_all = fwd_bed.intersect(genes_total_expressed, wa=True, nonamecheck=True)
    gene_reads_rev_all = rev_bed.intersect(genes_total_expressed, wa=True, nonamecheck=True)

    fwd_gene = get_IDs_from_bed(gene_reads_fwd_all)
    rev_gene = get_IDs_from_bed(gene_reads_rev_all)

    merged_gene = pd.concat([fwd_gene, rev_gene])
    merged_gene["ID"].str.replace('/1', '').str.replace('/2', '').drop_duplicates().to_csv(f"{aln_dir}/gene_reads.lst", header=None, index=False)
    print(colored("Done!", "green", attrs=['bold']))

    clock = time()
    print(f"{clock}\tChimeric reads pairs identification")
    with open(str(aln_dir + '/TE_reads.lst')) as f:
        list_TE_reads = f.read().splitlines()

    gene_reads = pd.read_csv(str(f"{aln_dir}/gene_reads.lst"), usecols=[0], names=['ID'])
    chim_reads = gene_reads[gene_reads['ID'].isin(list_TE_reads)].to_csv(f"{aln_dir}/chim_reads.lst", header=None, index=False)

    with open(str(f"{aln_dir}/chim_reads.lst")) as f:
        chim_list = f.read().splitlines()

    def overlap_reads(bed, list, output=None):
        intersect_TE = str(bed)
        intersect_TE = StringIO(intersect_TE)
        df = pd.read_table(intersect_TE, header=None, sep="\t", usecols=[0,1,2,3,4,5],names=['scaf', 'start', 'end', 'ID', 'dot', 'strand'])
        mod = df["ID"].str.replace('/1', '').str.replace('/2', '')
        df['ID_mod'] = mod
        if output is None:
            df[df['ID_mod'].isin(list)].to_csv(header=None, index=False, sep="\t", columns = ['scaf', 'start', 'end', 'ID', 'dot', 'strand'])
        else:
            df[df['ID_mod'].isin(list)].to_csv(str(aln_dir + '/' + output + '.bed'),header=None, index=False, sep="\t", columns = ['scaf', 'start', 'end', 'ID', 'dot', 'strand'])

    TE_reads_fwd = str('TE_reads_fwd')
    overlap_reads(TE_reads_fwd_all, chim_list, TE_reads_fwd)

    TE_reads_rev = str('TE_reads_rev')
    overlap_reads(TE_reads_rev_all, chim_list, TE_reads_rev)

    gene_reads_fwd = str('gene_reads_fwd')
    overlap_reads(gene_reads_fwd_all, chim_list, gene_reads_fwd)

    gene_reads_rev = str('gene_reads_rev')
    overlap_reads(gene_reads_rev_all, chim_list, gene_reads_rev)

    #Exons with chimeric reads
    accepted_hits = pybedtools.BedTool(str(f"{aln_dir}/accepted_hits.bed"))

    chim_reads_coord = str('chim_reads_coord')
    chim_reads = overlap_reads(accepted_hits, chim_list, chim_reads_coord)
    chim_reads = pybedtools.BedTool(str(f"{aln_dir}/chim_reads_coord.bed"))

    exon_file = pybedtools.BedTool(str(f"{tmp}/exon_file.bed"))
    chimeric_exons = exon_file.intersect(chim_reads, wa=True, nonamecheck=True)
    chimeric_exons = dropdup_bed(chimeric_exons)
    with open(str(f"{aln_dir}/chim_exons.bed"), 'w') as f:
        print(str(chimeric_exons), file = f)

    chimeric_TEs = expressed_TEs.intersect(chim_reads, wa=True, nonamecheck=True, F=str(args.overlap))
    chimeric_TEs = dropdup_bed(chimeric_TEs)
    with open(str(f"{aln_dir}/chimeric_TEs.bed"), 'w') as f:
        print(str(chimeric_TEs), file = f)

    print(colored("Done!", "green", attrs=['bold']))
