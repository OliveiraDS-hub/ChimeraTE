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
parser.parse_args()
args = parser.parse_args()

print(f"/==================== Project {args.project} ====================\\")

print(f"\nTranscripts file ==>\t{args.transcripts}")
print(f"TE file ==>\t{args.te}")
if args.assembly:
    print(f"Transcriptome assembly ==> ON\n")
else:
    print(f"Transcriptome assembly ==> OFF\n")

out_genome = str(args.transcripts).replace(".fasta","").replace(".fas","").replace(".fa","")
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

def mate_spec_IDs(bedfile, mate):
    df = pd.read_csv(str(bedfile), header=None, sep="\t", usecols=[3],names=['read_ID'])
    df = df.to_csv(header=None, index=False, sep="\t"); df = StringIO(df)
    pattern = str('/' + str(mate))
    reads = ''
    for line in df:
        if re.search(pattern, line):
            reads += line.replace(str(pattern), '')
    return reads

def import_csv(dataframe):
    df = pd.read_csv(str(dataframe), header=None, sep="\t", usecols=[6,7,8,9,10,11],names=['scaf', 'start', 'end', 'ID', 'dot', 'strand'])
    return df

def overlap(bed, list, out):
    df = pd.read_csv(str(bed), header=None, sep="\t", usecols=[0,1,2,3,4,5],names=['scaf', 'start', 'end', 'ID', 'dot', 'strand'])
    df[['ID']] = df[['ID']].replace('/1', '', regex=True).replace('/2', '', regex=True)
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

def rm_file_or_dir(f):
    if os.path.exists(f):
        if os.path.isfile(f):
            os.remove(f)
        elif os.path.isdir(f):
            if not os.listdir(f):
                os.rmdir(f)

def check_file(f):
    if os.path.exists(f):
        if os.stat(str(f)).st_size > 0:
            return True
        else:
            return False
    else:
        return False

out_dir = create_dir(str(f"{mydir}/projects/{args.project}"))
tmp = create_dir(str(f"{mydir}/projects/{args.project}/tmp"))
create_dir(str(f"{out_dir}/index"))

###import a script
import mode2_prep

###import a function from a script
from mode2_alignment import alignment_func
from mode2_chim_transcripts import multicore_chimeras
from mode2_chim_transcripts import prep_data
from mode2_chim_transcripts import merging_transc
from mode2_chim_transcripts import expression
from mode2_assembly import transcriptome_assembly
from mode2_assembly import singleton_crossing
from mode2_replicability import chim_reads_rep
from mode2_replicability import trasnc_rep

for index, row in input.iterrows():
    mate1 = row['mate1']
    mate2 = row['mate2']
    group = row['group']
    out_group = create_dir(str(f"{out_dir}/{group}"))
    global aln_dir
    print(f"Running analysis with ------------------------------------------> {group}\n")
    aln_dir = create_dir(str(f"{out_dir}/{group}/alignment"))
    create_dir(str(f"{out_dir}/{group}/alignment/fpkm_counts"))

    ##Perform bowtie2 alignment and identify chimeric reads
    alignment_func(out_dir, aln_dir, mate1, mate2)

    ##Perform STAR alignment and identify chimeric reads
    prep_data()
    multicore_chimeras()

    if check_file(str(f"{out_group}/chimTEs_raw.tsv")) == True:
        merging_transc()
    if check_file(str(f"{out_group}/chimTEs_final.tsv")) == True:
        expression()

    ###Assembly
    if args.assembly:
        trinity_out = create_dir(str(f"{out_dir}/{group}/trinity_out"))
        transcriptome_assembly()
        singleton_crossing()

chim_reads_rep()
if args.assembly:
    trasnc_rep()
print(colored("ChimeraTE has finished successfullly!", "green", attrs=['bold']))
