import sys
sys.dont_write_bytecode = True
import argparse
from argparse import ArgumentParser,SUPPRESS
from termcolor import colored
from io import StringIO
import pandas as pd
import subprocess
import pybedtools
import contextlib
import textwrap
import datetime
import shutil
import glob
import csv
import os
import re

from scripts.mode2_prep import *
from scripts.mode2_alignment import *
from scripts.mode2_chim_transcripts import *
from scripts.mode2_assembly import *
from scripts.mode2_replicability import *

print(colored("   ________    _                         ","white")+ colored(' ____________', 'red', attrs=['bold']))
print(colored("  / ____/ /_  (_)___ ___  ___  _________","white")+ colored(' /_  __/ ____/', 'red', attrs=['bold']))
print(colored(" / /   / __ \/ / __ ` __\/ _ \/ ___/ __ `","white") + colored(' / / / __/', 'red', attrs=['bold']))
print(colored("/ /___/ / / / / / / / / /  __/ /  / /_/ /","white") + colored('/ / / /___', 'red', attrs=['bold']))
print(colored("\____/_/ /_/_/_/ /_/ /_/\___/_/   \__,_/", "white") + colored('/_/ /_____/', 'red', attrs=['bold']))
print(colored("-. .-.   .-. .-.   .-. .-.   .-. .-.   .", "white") + colored('-. .-.   .-. .', 'red', attrs=['bold']))
print(colored("||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|", "white") + colored('||\|||\ /|||\|', 'red', attrs=['bold']))
print(colored("|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\||", "white") + colored('|/ \|||\|||/ \ ', 'red', attrs=['bold']))
print(colored('~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-', "white") + colored('`   `-~ `-`   ', 'red', attrs=['bold']))
print("Version 1.2")

sys.path.insert(1, 'scripts/')
parser = argparse.ArgumentParser(description='ChimeraTE Mode 2: The genome-blinded approach to detect chimeric transcripts with RNA-seq data.', usage=SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter, epilog=textwrap.dedent('''Citation: Oliveira, D. S., Fablet, M., Larue, A., Vallier, A., Carareto, C. M., Rebollo, R., & Vieira, C. (2023). ChimeraTE: a pipeline to detect chimeric transcripts derived from genes and transposable elements. Nucleic Acids Research, 51(18), 9764-9784.'''))
parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')
required.add_argument('--input', help='Paired-end files and their respective group/replicate', required=True, type=str, metavar = "")
required.add_argument('--project', help='Directory name with output data', required=True, type=str, metavar = "")
required.add_argument('--te', help='Fasta file containing TE information', required=True, type=str, metavar = "")
required.add_argument('--transcripts', help='Fasta file containing gene information', required=True, type=str, metavar = "")
required.add_argument('--strand', choices=['rf-stranded','fwd-stranded'], required=True, help='Define the strandness direction of the RNA-seq. Two options: \"rf-stranded\" OR \"fwd-stranded\"', type=str, metavar = "")

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('--coverage', help='Minimum coverage (mean between replicates default 2 for chimeric transcripts detection)', required=False, type=str, default=2, metavar = "")
optional.add_argument('--fpkm', help='Minimum fpkm to consider a gene as expressed (default 1)', required=False, type=str, default=1, metavar = "")
optional.add_argument('--replicate', help='Minimum recurrency of chimeric transcripts between RNA-seq replicates (default 2)', required=False, type=str, default=2, metavar = "")
optional.add_argument('--threads', help='Number of threads (default 6)', required=False, type=str, default=6, metavar = "")
optional.add_argument('--assembly', help='Search for chimeric transcript with transcriptome assembly with Trinity', required=False, action = 'store_true')
optional.add_argument('--ref_TEs', help='"species" database used by RepeatMasker (flies, human, mouse, arabidopsis; or a built TE library in fasta format)', required=False, type=str, default=str('flies'), metavar = "")
optional.add_argument('--ram', help='Minimum RAM memory in Gbytes (default 8)', required=False, type=str, default=8, metavar = "")
optional.add_argument('--overlap', help='Minimum overlap between chimeric reads and TE insertions (default 0.50)', required=False, type=float, default=0.50, metavar = "")
optional.add_argument('--TE_length', help='Minimum TE length to keep it from RepeatMasker output (default = 80bp)', required=False, type=int, default=80, metavar = "")
optional.add_argument('--identity', help='Minimum identity between de novo assembled transcripts and reference transcripts (default = 80)', required=False, type=int, default=80, metavar = "")
optional.add_argument('--index', help='Folder with bowtie2 index from transcriptome', required=False, type=str, metavar = "")

parser.parse_args()
args = parser.parse_args()

print(f"/==================== Project {args.project} ====================\\")

print(f"\nTranscripts file ==>\t{args.transcripts}")
print(f"TE file ==>\t{args.te}")
if args.assembly:
    print(f"Transcriptome assembly ==> ON\n")
else:
    print(f"Transcriptome assembly ==> OFF\n")

from util.parser_mode2 import *

out_genome = str(args.transcripts).replace(".fasta","").replace(".fas","").replace(".fa","")
input = pd.read_csv(args.input, header=None, sep="\t", usecols=[0,1,2],names=['mate1', 'mate2', 'group'])
mydir = str(os.getcwd())

out_dir = create_dir(str(f"{mydir}/projects/{args.project}"))
create_dir(f'{mydir}/index')
tmp = create_dir(str(f"{mydir}/projects/{args.project}/tmp"))
create_dir(str(f"{out_dir}/index"))

### Create bowtie2 indexes
bowtie_indexer(args.te, args.transcripts, args.threads)

for index, row in input.iterrows():
    mate1 = row['mate1']
    mate2 = row['mate2']
    group = row['group']
    out_group = create_dir(f"{out_dir}/{group}")

    print(f"Running analysis with ------------------------------------------> {group}\n")
    aln_dir = create_dir(f"{out_dir}/{group}/alignment")
    create_dir(f"{out_dir}/{group}/alignment/fpkm_counts")

    ##Perform bowtie2 alignment and identify chimeric reads
    alignment_func(out_dir, args.index, aln_dir, mate1, mate2, args.threads, args.strand, args.transcripts)

    ##Perform STAR alignment and identify chimeric reads
    prep_data()
    multicore_chimeras(args.threads)

    if check_file(f"{out_group}/chimTEs_raw.tsv") == True:
        merging_transc()
    if check_file(f"{out_group}/chimTEs_final.tsv") == True:
        expression()

    ###Assembly
    if args.assembly:
        trinity_out = create_dir(f"{out_dir}/{group}/trinity_out")
        transcriptome_assembly(mate1, mate2, args.strand, args.ram, args.threads, args.ref_TEs, args.TE_length, args.overlap)
        singleton_crossing(args.transcripts, args.identity)

chim_reads_rep()
if args.assembly:
    trasnc_rep()
print(colored("ChimeraTE has finished successfullly!", "green", attrs=['bold']))
