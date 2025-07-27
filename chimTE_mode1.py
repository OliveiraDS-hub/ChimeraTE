import sys ;sys.dont_write_bytecode = True
from argparse import ArgumentParser,SUPPRESS
from termcolor import colored
from io import StringIO
import pyranges as pr
import pandas as pd
import subprocess
import contextlib
import argparse
import datetime
import textwrap
import shutil
import glob
import csv
import os
import re

from util.parser import *
from util.intersector import *
from scripts.mode1_prep_data import *
from scripts.mode1_alignment import *
from scripts.mode1_te_initiated import *
from scripts.mode1_te_terminated import *
from scripts.mode1_te_exonized import *
from scripts.mode1_replicability import *

print(colored("   ________    _                         ","white")+ colored(' ____________', 'red', attrs=['bold']))
print(colored("  / ____/ /_  (_)___ ___  ___  _________","white")+ colored(' /_  __/ ____/', 'red', attrs=['bold']))
print(colored(" / /   / __ \/ / __ ` __\/ _ \/ ___/ __ `","white") + colored(' / / / __/', 'red', attrs=['bold']))
print(colored("/ /___/ / / / / / / / / /  __/ /  / /_/ /","white") + colored('/ / / /___', 'red', attrs=['bold']))
print(colored("\____/_/ /_/_/_/ /_/ /_/\___/_/   \__,_/", "white") + colored('/_/ /_____/', 'red', attrs=['bold']))
print(colored("-. .-.   .-. .-.   .-. .-.   .-. .-.   .", "white") + colored('-. .-.   .-. .', 'red', attrs=['bold']))
print(colored("||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|", "white") + colored('||\|||\ /|||\|', 'red', attrs=['bold']))
print(colored("|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\||", "white") + colored('|/ \|||\|||/ \ ', 'red', attrs=['bold']))
print(colored('~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-', "white") + colored('`   `-~ `-`   ', 'red', attrs=['bold']))
print("Version 2.0.0\n")

sys.path.insert(1, 'scripts/')
parser = argparse.ArgumentParser(description='ChimeraTE Mode 1: The genome-guided approach to detect chimeric transcripts with RNA-seq data.', usage=SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter, epilog=textwrap.dedent('''Citation: Oliveira, D. S., Fablet, M., Larue, A., Vallier, A., Carareto, C. M., Rebollo, R., & Vieira, C. (2023). ChimeraTE: a pipeline to detect chimeric transcripts derived from genes and transposable elements. Nucleic Acids Research, 51(18), 9764-9784.'''))
parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')
required.add_argument('--genome', help='Genome in fasta', required=True, type=str, metavar = "")
required.add_argument('--input', help='Paired-end files and their respective group/replicate', required=True, type=str, metavar = "")
required.add_argument('--project', help='Directory name with output data', required=True, type=str, metavar = "")
required.add_argument('--te', help='GTF file containing TE information', required=True, type=str, metavar = "")
required.add_argument('--gene', help='GTF file containing gene information', required=True, type=str, metavar = "")
required.add_argument('--strand', choices=['rf-stranded','fwd-stranded'], required=True, help='Define the strandness direction of the RNA-seq. Two options: \"rf-stranded\" OR \"fwd-stranded\"', type=str, metavar = "")

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('--chimera', choices=['TE-initiated','TE-exonized', 'TE-terminated'], help='Identify specific type of chimera: \"TE-initiated\" OR \"TE-exonized\" OR \"TE-terminated\"', required=False, type=str, metavar = "")
optional.add_argument('--window', help='Upstream and downstream window size (default = 3000)', required=False, type=str, default=3000, metavar = "")
optional.add_argument('--replicate', help='Minimum recurrency of chimeric transcripts between RNA-seq replicates (default = 2)', required=False, type=str, default=2, metavar = "")
optional.add_argument('--coverage', help='Minimum coverage as the mean between replicates for chimeric transcripts detection', required=False, type=str, default=2, metavar = "")
optional.add_argument('--fpkm', help='Minimum fpkm to consider a gene as expressed (default = 1)', required=False, type=str, default=1, metavar = "")
optional.add_argument('--threads', help='Number of threads (default = 6)', required=False, type=int, default=6, metavar = "")
optional.add_argument('--overlap', help='Minimum overlap between chimeric reads and TE insertions (default = 0.50)', required=False, type=float, default=0.50, metavar = "")
optional.add_argument('--index', help='Absolute path to STAR index', required=False, type=str, metavar = "")
parser.parse_args()
args = parser.parse_args()
print(f"/==================== Project {args.project} ====================\\")

out_genome = str(args.genome).replace(".fasta","").replace(".fas","").replace(".fa","")
input = pd.read_csv(args.input, header=None, sep="\t", usecols=[0,1,2],names=['mate1', 'mate2', 'group'])
mydir = str(os.getcwd())
print(out_genome)
out_dir = create_dir(str(f"{mydir}/projects/{args.project}"))
tmp = create_dir(str(f"{mydir}/projects/{args.project}/tmp"))


### Check GTF annotations
clock = time()
print(f"{clock}\tChecking gene and TE annotations")

features_GTF = ['gene', 'exon']
print(f"GTF GENE\n{args.gene} contains:", end = '\n')
for feat in features_GTF:
    res_test = test_GTF_feature(str(tmp), str(args.gene), str(feat))
    if res_test == False:
        print(str(f"ERROR: Bad GTF format\n{args.gene.name} does not contain coordinates of {feat}s! The 3rd column must contain \"{feat}\"\tExiting..."))
        exit()

## Check gene IDs
n_genes = gene_IDs(tmp)
print(f"First two gene IDs: {n_genes}")

print(f"\nGTF TE\n{args.te} contains:", end = '\n')
count_TE_families(str(args.te))



### Create STAR index and bed files
annotation_manager(tmp, args.index, out_dir, args.genome, args.te, args.threads)

for index, row in input.iterrows():
    global aln_dir
    mate1 = row['mate1']
    mate2 = row['mate2']
    group = row['group']

    out_group = create_dir(f"{out_dir}/{group}")
    aln_dir = create_dir(f"{out_dir}/{group}/alignment")

    ###Perform STAR alignment and identify chimeric reads
    aln_check = alignment_func(out_dir,group,aln_dir,mate1,mate2, args.threads, args.strand, tmp, args.fpkm)
    # aln_check = True
    if aln_check:
        if args.chimera is None or args.chimera == "TE-initiated":
            ###Search for TE-initiated transcripts
            if not os.path.exists(f"{out_group}/TE-initiated-{group}.tsv"):
                itsct_TEs_up = te_init(aln_dir, args.window)
                multicore_process_init(itsct_TEs_up, args.threads)
            else:
                print(f"{out_group}/TE-initiated-{group}.tsv has been found!\tSkipping the identification of TE-initiated transcripts...")

        if args.chimera is None or args.chimera == "TE-terminated":
            ###Search for TE-terminated transcripts
            if not os.path.exists(f"{out_group}/TE-terminated-{group}.tsv"):
                itsct_TEs_down = te_term(aln_dir, args.window)
                multicore_process_term(itsct_TEs_down, args.threads)
            else:
                print(f"{out_group}/TE-terminated-{group}.tsv has been found!\tSkipping the identification of TE-initiated transcripts...")

        if args.chimera is None or args.chimera == "TE-exonized":
            ### Search for TE-exonized transcripts - embedded
            te_emb_ids, te_emb_mbed, embedded_exons = prep_embedded(aln_dir)
            multicore_emb_exon(te_emb_ids, te_emb_mbed, args.threads, args.overlap)

            ### Search for TE-exonized transcripts - overlapped
            overlap_TEs_exons_mbed, geneIDs_overl = prep_overlapped(aln_dir, te_emb_mbed)
            multicore_overl_exon(geneIDs_overl, overlap_TEs_exons_mbed, args.threads, args.overlap)

            ###Search for TE-exonized transcripts - intronic
            geneids_intron, exons_genes_TEs_intron, genes_TE_intron_mbed = prep_intronic(aln_dir)
            multicore_intron_exon(geneids_intron, exons_genes_TEs_intron, genes_TE_intron_mbed, args.threads , args.overlap)

        # if os.path.exists(f"{out_group}/TE-exonized-{group}.tsv") == True:
        #     copy(f"{out_group}/TE-exonized-{group}.tsv", tmp)
        # else:
        #     print(f"Skipping the identification of TE-exonized transcripts...")

###Checking for replicability of chimeric transcripts in RNA-seq samples
replicability(args.replicate, args.coverage)

print(colored("ChimeraTE has finished successfullly!", "green", attrs=['bold']))