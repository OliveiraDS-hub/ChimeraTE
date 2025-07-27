import subprocess
import pybedtools
import pandas as pd
import datetime
import shutil
import os
from io import StringIO
import glob

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

def dropdup_bed(bedtool_obj):
    # Dump to a temp file and read its content
    temp_file = bedtool_obj.saveas()
    with open(temp_file.fn) as f:
        bed_str = f.read()

    if not bed_str.strip():
        return pybedtools.BedTool('')  # Return empty BedTool

    df = pd.read_csv(StringIO(bed_str), sep="\t", header=None).drop_duplicates()
    return pybedtools.BedTool.from_dataframe(df)


def get_IDs_from_bed(bed_file):
    if isinstance(bed_file, pd.DataFrame):
        # Only select the 4th column
        ID_column = bed_file.iloc[:, [3]].copy()
        ID_column.columns = ["Name"]
    elif os.path.exists(bed_file):
        # Read only the 4th column from file
        ID_column = pd.read_csv(
            bed_file,
            sep="\t",
            header=None,
            usecols=[3],         # only load column 4
            names=["Name"]       # assign column name directly
        )
    return ID_column

def import_csv(dataframe):
    df = pd.read_csv(str(dataframe), header=None, sep="\t", usecols=[6,7,8,9,10,11],names=['scaf', 'start', 'end', 'ID', 'dot', 'strand'])
    return df

# def overlap(bed, list, out):
#     df = pd.read_csv(str(bed), header=None, sep="\t", usecols=[0,1,2,3,4,5],names=['scaf', 'start', 'end', 'ID', 'dot', 'strand'])
#     df[df['ID'].isin(list)].to_csv(str(out), header=None, index=False, sep="\t")

def overlap(bed, ids, out):
    id_set = set(ids)
    with open(bed, 'r') as infile, open(out, 'w') as outfile:
        for line in infile:
            if line.split('\t')[3] in id_set:
                outfile.write(line)

def time():
    x = datetime.datetime.now()
    clock = str("[" + str(x.strftime("%A")) + " " + str(x.day) + "/" + str(x.month) +  "/" + str(x.year) + " - " + str(x.strftime("%H")) + "h:" + str(x.strftime("%M")) + "]")
    return clock

def copy(file, folder):
    if os.path.exists(file):
        shutil.copy(file, folder)

def remove(dir):
    list = ["rev*", "fwd*", "accepted_hits.bed", "*out.bam", "gene_reads.lst"]
    for pattern in list:
        for file in glob.glob(str(f"{dir}/{pattern}")):
            os.remove(file)