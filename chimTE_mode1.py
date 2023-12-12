import sys
sys.dont_write_bytecode = True
import argparse
import datetime
from termcolor import colored
import shutil
import subprocess
import glob, os
import pandas as pd
import re
import csv
import textwrap
from argparse import ArgumentParser,SUPPRESS
from io import StringIO
import pybedtools
import contextlib
from termcolor import colored
import sys

print(colored("   ________    _                         ","white")+ colored(' ____________', 'red', attrs=['bold']))
print(colored("  / ____/ /_  (_)___ ___  ___  _________","white")+ colored(' /_  __/ ____/', 'red', attrs=['bold']))
print(colored(" / /   / __ \/ / __ ` __\/ _ \/ ___/ __ `","white") + colored(' / / / __/', 'red', attrs=['bold']))
print(colored("/ /___/ / / / / / / / / /  __/ /  / /_/ /","white") + colored('/ / / /___', 'red', attrs=['bold']))
print(colored("\____/_/ /_/_/_/ /_/ /_/\___/_/   \__,_/", "white") + colored('/_/ /_____/', 'red', attrs=['bold']))
print(colored("-. .-.   .-. .-.   .-. .-.   .-. .-.   .", "white") + colored('-. .-.   .-. .', 'red', attrs=['bold']))
print(colored("||\|||\ /|||\|||\ /|||\|||\ /|||\|||\ /|", "white") + colored('||\|||\ /|||\|', 'red', attrs=['bold']))
print(colored("|/ \|||\|||/ \|||\|||/ \|||\|||/ \|||\||", "white") + colored('|/ \|||\|||/ \ ', 'red', attrs=['bold']))
print(colored('~   `-~ `-`   `-~ `-`   `-~ `-~   `-~ `-', "white") + colored('`   `-~ `-`   ', 'red', attrs=['bold']))
print("Version 1.1.1")

sys.path.insert(1, 'scripts/')
parser = argparse.ArgumentParser(description='ChimeraTE Mode 1: The genome-guided approach to detect chimeric transcripts with RNA-seq data.', usage=SUPPRESS, formatter_class=argparse.RawDescriptionHelpFormatter, epilog=textwrap.dedent('''Citation: Oliveira, D. S., et al. (2022). ChimeraTE: A pipeline to detect chimeric transcripts derived from genes and transposable elements. bioRxiv, 2022-09.'''))
parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')
required.add_argument('--genome', help='Genome in fasta', required=True, type=str, metavar = "")
required.add_argument('--input', help='Paired-end files and their respective group/replicate', required=True, type=str, metavar = "")
required.add_argument('--project', help='Directory name with output data', required=True, type=str, metavar = "")
required.add_argument('--te', help='GTF file containing TE information', required=True, type=str, metavar = "")
required.add_argument('--gene', help='GTF file containing gene information', required=True, type=argparse.FileType('r'), metavar = "")
required.add_argument('--strand', choices=['rf-stranded','fwd-stranded'], required=True, help='Define the strandness direction of the RNA-seq. Two options: \"rf-stranded\" OR \"fwd-stranded\"', type=str, metavar = "")

optional = parser.add_argument_group('Optional arguments')
optional.add_argument('--window', help='Upstream and downstream window size (default = 3000)', required=False, type=str, default=3000, metavar = "")
optional.add_argument('--replicate', help='Minimum recurrency of chimeric transcripts between RNA-seq replicates (default 2)', required=False, type=str, default=2, metavar = "")
optional.add_argument('--coverage', help='Minimum coverage (mean between replicates default 2 for chimeric transcripts detection)', required=False, type=str, default=2, metavar = "")
optional.add_argument('--fpkm', help='Minimum fpkm to consider a gene as expressed (default 1)', required=False, type=str, default=1, metavar = "")
optional.add_argument('--threads', help='Number of threads (default 6', required=False, type=str, default=6, metavar = "")
optional.add_argument('--overlap', help='Minimum overlap between chimeric reads and TE insertions (default 0.50)', required=False, type=float, default=0.50, metavar = "")
optional.add_argument('--index', help='Absolute path to STAR index', required=False, type=str, metavar = "")
parser.parse_args()
args = parser.parse_args()

out_genome = str(args.genome).replace(".fasta","").replace(".fas","").replace(".fa","")
input = pd.read_csv(args.input, header=None, sep="\t", usecols=[0,1,2],names=['mate1', 'mate2', 'group'])
mydir = str(os.getcwd())

def create_dir(path):
    if not os.path.isdir(path):
        os.mkdir(path)
    return path

def samt_index(bam_file):
    subprocess.call(['samtools', 'index', str(bam_file)])

def intersection(file1, file2, output=None, prop=None):
    if prop is not None:
        with open(str(output), 'w') as tmp_var:
            subprocess.call(['bedtools', 'intersect', '-a', str(file1), '-b', str(file2), '-wa', '-wb', '-f', str(prop), '-nonamecheck'], stdout=tmp_var)
        tmp_var.close
    else:
        with open(str(output), 'w') as tmp_var:
            subprocess.call(['bedtools', 'intersect', '-a', str(file1), '-b', str(file2), '-wa', '-wb', '-nonamecheck'], stdout=tmp_var)
        tmp_var.close

def pybedtools_intersection(file1, file2, prop = None):
    bed1 = pybedtools.BedTool(str(file1))
    if type(file2) == str:
        bed2 = pybedtools.BedTool(str(file2))
    else:
        bed2 = file2
    if prop is not None:
        intersect = bed1.intersect(bed2, wa=True, nonamecheck=True, f=float(prop))
    else:
        intersect = bed1.intersect(bed2, wa=True, nonamecheck=True)
    return intersect

def dropdup_bed(intersection_bed):
    bed2string = str(intersection_bed)
    bed2string = StringIO(bed2string)
    pd_df = pd.read_table(bed2string, header=None, sep="\t").drop_duplicates().to_csv(header=None, index=False, sep="\t")
    return pd_df

def get_IDs_from_bed(bed_file):
    bed2string = str(bed_file)
    bed2string = StringIO(bed2string)
    ID_column = pd.read_csv(bed2string, sep="\t", usecols=[3],names=['ID'])
    return ID_column

def import_csv(dataframe):
    df = pd.read_csv(str(dataframe), header=None, sep="\t", usecols=[6,7,8,9,10,11],names=['scaf', 'start', 'end', 'ID', 'dot', 'strand'])
    return df

def overlap(bed, list, out):
    df = pd.read_csv(str(bed), header=None, sep="\t", usecols=[0,1,2,3,4,5],names=['scaf', 'start', 'end', 'ID', 'dot', 'strand'])
    df[df['ID'].isin(list)].to_csv(str(out), header=None, index=False, sep="\t")

def time():
    x = datetime.datetime.now()
    clock = str("[" + str(x.strftime("%A")) + " " + str(x.day) + "/" + str(x.month) +  "/" + str(x.year) + " - " + str(x.strftime("%H")) + "h:" + str(x.strftime("%M")) + "]")
    return clock

def copy(file, folder):
    if os.path.exists(file):
        shutil.copy(file, folder)

def remove():
    list = ["rev*", "fwd*", "accepted_hits.bed", "*out.bam", "gene_reads.lst"]
    for pattern in list:
        for file in glob.glob(str(aln_dir) + '/' + str(pattern)):
            os.remove(file)

out_dir = create_dir(str(mydir+'/projects/'+str(args.project)))
tmp = create_dir(str(mydir+'/projects/'+str(args.project)+'/tmp'))

with open(str(tmp + '/gtf_file.gtf'), 'w', encoding='utf-8') as f:
    for line in (args.gene):
        if not (line.startswith('#')):
            f.write(line)
f.close()

## Create STAR index and bed files
import mode1_prep_data
from mode1_alignment import alignment_func
from mode1_te_initiated import te_init
from mode1_te_initiated import init_chimeras
from mode1_te_initiated import multicore_process_init
from mode1_te_terminated import te_term
from mode1_te_terminated import multicore_process_term
from mode1_te_exonized import te_exon_embedded
from mode1_te_exonized import prep_overlapped
from mode1_te_exonized import prep_intronic
from mode1_te_exonized import multicore_process_exon

for index, row in input.iterrows():
    mate1 = row['mate1']
    mate2 = row['mate2']
    group = row['group']
    out_group = create_dir(str(out_dir + '/' + group))
    global aln_dir
    aln_dir = create_dir(str(out_dir + '/' + group + '/alignment'))

    ###Perform STAR alignment and identify chimeric reads
    alignment_func(out_dir,group,aln_dir,mate1,mate2)

    ###Search for TE-initiated transcripts
    te_init(aln_dir, group, out_group)
    multicore_process_init()
    init_output = str(out_group + "/TE-initiated-" + str(group) + ".tsv")
    if os.path.exists(init_output) == True:
        copy(str(init_output), tmp)

    ###Search for TE-terminated transcripts
    te_term(aln_dir, group, out_group)
    multicore_process_term()
    term_output = str(out_group + "/TE-terminated-" + str(group) + ".tsv")
    if os.path.exists(term_output) == True:
        copy(str(term_output), tmp)

    ##Search for TE-exonized transcripts
    te_exon_embedded(aln_dir, group, out_group)
    prep_overlapped(aln_dir, group, out_group)
    prep_intronic(aln_dir, group, out_group)
    multicore_process_exon()
    exon_output = str(out_group + "/TE-exonized-" + str(group) + ".tsv")
    if os.path.exists(exon_output) == True:
        copy(str(exon_output), tmp)

    ##Removing tmp files
    remove()

#Checking for replicability of chimeric transcripts in RNA-seq samples
import mode1_replicability









#
