import os
import subprocess
import pybedtools
from io import StringIO
import pandas as pd
import datetime
import shutil
import glob
import re


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
