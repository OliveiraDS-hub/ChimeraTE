import pyranges as pr
import pandas as pd
import numpy as np
import subprocess
import argparse
import re
import os

parser = argparse.ArgumentParser(description='')
parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')
required.add_argument('--input', help='Paired-end files and their respective group/replicate', required=True, type=str, metavar = "")
required.add_argument('--index', help='Directory STAR index', required=True, type=str, metavar = "")
required.add_argument('--output', help='prefix for output files', required=True, type=str, metavar = "")
required.add_argument('--gff', help='gff with gene annotation', required=True, type=str, metavar = "")
parser.parse_args()
args = parser.parse_args()


def getTPM(df, out_prefix):
    stringtie_df = pd.read_csv(df, header=None, 
                               sep="\t",
                               usecols=[0,2,3,4,6,7,8],
                               names=['scaf', 'feature', 'start','end', 'strand', 'dot','ID'],
                               comment='#')
    stringtie_df['transcriptID'] = stringtie_df['ID'].str.split(';').str[1]
    stringtie_df['transcriptID'] = stringtie_df['transcriptID'].str.replace("\"", '').str.replace("transcript_id", '').str.replace(" ", '')
    stringtie_transcripts = stringtie_df[stringtie_df["feature"] == "transcript"]

    stringtie_transcripts['TPM'] = stringtie_transcripts['ID'].str.extract(r'TPM "([^"]+)"').astype(float)
    stringtie_transcripts['ref_gene'] = stringtie_transcripts['ID'].str.extract(r'ref_gene_name "([^"]+)"').astype(str)

    # Convert string "nan" to real NaN
    stringtie_transcripts['ref_gene'].replace("nan", pd.NA, inplace=True)

    # Extract gene_id from transcriptID (STRG.22.1 → STRG.22)
    stringtie_transcripts['gene_id'] = stringtie_transcripts['transcriptID'].str.extract(r'(STRG\.\d+)')

    # Fill NaNs in ref_gene per gene_id
    stringtie_transcripts['ref_gene'] = (
        stringtie_transcripts.groupby('gene_id')['ref_gene']
        .transform(lambda x: x.fillna(method='ffill').fillna(method='bfill'))
    )

    # Keep only ref_gene, transcriptID, TPM
    tpm_transcripts = stringtie_transcripts[['ref_gene', 'transcriptID', 'TPM']]
    tpm_transcripts.to_csv(f"{out_prefix}_stringtie_TPM_values.csv", index=False)


    ## extract bed exons
    stringtie_exons = stringtie_df[stringtie_df["feature"] == "exon"]
    stringtie_exons['TPM'] = stringtie_exons['ID'].str.extract(r'TPM "([^"]+)"').astype(float)
    stringtie_exons['transcriptID'] = stringtie_exons['ID'].str.split(';').str[1]
    stringtie_exons['transcriptID'] = stringtie_exons['transcriptID'].str.replace("\"", '').str.replace("transcript_id", '').str.replace(" ", '')
    bed_exons = stringtie_exons[['scaf', 'start', 'end', 'transcriptID', 'dot', 'strand']]
    return bed_exons

def import_chimeras(df):
    chim_df = pd.read_csv(df, sep="\t",
                          usecols=[0,3,5])
    return chim_df

def overlap_fraction_intersection(a_bed_file, b_bed_file, fraction=1.0):
    """
    Intersect A and B using pyranges, keeping only A intervals where the fraction
    of A overlapped by B is >= `fraction`.

    Equivalent to: bedtools intersect -a A -b B -wa -wb -f <fraction>
    """
    # Read BED files
    if isinstance(a_bed_file, pd.DataFrame):
        if len(a_bed_file.columns) == 3:
            a_bed_file.columns = ["Chromosome", "Start", "End"]
            print(a_bed_file)
            a_bed_file.to_csv("TEs.bed", sep = "\t", index = False, header=None)
        a = pr.PyRanges(a_bed_file)
    else:
        a = pr.read_bed(a_bed_file)
    
    if isinstance(b_bed_file, pd.DataFrame):
        if len(b_bed_file.columns) == 6:
            b_bed_file.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
            b_bed_file.to_csv("exons.bed", sep = "\t", index = False, header=None)
        b = pr.PyRanges(b_bed_file)
    else:
        b = pr.read_bed(b_bed_file)

    # Join and compute overlap
    joined = a.join(b, suffix="_b", how="left", report_overlap=True)

    df = joined.df.copy()
    df["A_length"] = df["End"] - df["Start"]
    df["overlap_fraction"] = df["Overlap"] / df["A_length"]

    # Filter for required fraction
    df = df[df["overlap_fraction"] >= fraction]

    return df



#### Merge replicates
# Read the TSV file
df = pd.read_csv(args.input, sep="\t", header=None, names=["R1", "R2", "rep"])

# Get prefix from first R1 file
first_file = os.path.basename(df["R1"].iloc[0])
prefix = first_file.split("_P")[0]

## Check if merged files already exists
merged_R1_exists = os.path.exists("merged_R1.fastq.gz") or os.path.exists("merged_R1.fastq")
merged_R2_exists = os.path.exists("merged_R2.fastq.gz") or os.path.exists("merged_R2.fastq")

if not merged_R1_exists and not merged_R2_exists:
    # Get all R1 and R2 file paths
    R1_files = df.iloc[:, 0].tolist()  # First column = R1
    R2_files = df.iloc[:, 1].tolist()  # Second column = R2

    # Detect compression from first R1 file
    compressed = R1_files[0].endswith(".gz")
    output_R1 = "merged_R1.fastq.gz" if compressed else "merged_R1.fastq"
    output_R2 = "merged_R2.fastq.gz" if compressed else "merged_R2.fastq"

    # Merge R1 files
    print(f"Merging {len(R1_files)} R1 files into {output_R1}...")
    if compressed:
        subprocess.run(["cat"] + R1_files, stdout=open(output_R1, "wb"), check=True)
    else:
        subprocess.run(["cat"] + R1_files, stdout=open(output_R1, "w"), check=True)

    # Merge R2 files
    print(f"Merging {len(R2_files)} R2 files into {output_R2}...")
    if compressed:
        subprocess.run(["cat"] + R2_files, stdout=open(output_R2, "wb"), check=True)
    else:
        subprocess.run(["cat"] + R2_files, stdout=open(output_R2, "w"), check=True)
else:
    print(f"Merged files have been found!\tSkipping...")

### STAR alignment
if not os.path.exists(f"{args.output}_Aligned.sortedByCoord.out.bam"):
    subprocess.run(f'STAR --genomeDir {args.index} --runThreadN 16 --readFilesCommand zcat --readFilesIn merged_R1.fastq.gz merged_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {args.output}_ --outTmpDir STARtmp', shell=True)
else:
    print("Alignment bam has been found!")


### Stringtie expression
subprocess.run(f'stringtie {args.output}_Aligned.sortedByCoord.out.bam -G {args.gff} -o {args.output}_gene_sample.gtf -A {args.output}_gene_abundances.tsv -B -p 30', shell=True)

## Create bed for exon positions
bed_exons = getTPM(f"{args.output}_gene_sample.gtf", str(args.output))
## Import stringtie outputs
stringtie_TPM_df = pd.read_csv(f"{args.output}_stringtie_TPM_values.csv")
stringtie_abund = pd.read_csv(f"{args.output}_gene_abundances.tsv", usecols=[0,8], sep ="\t")

### Store the results
gene_ids = []
total_gene_TPMs = []
chimeric_TPMs = []
contributions = []

### Add contribution to chimTE's output
chimTE_out_df = pd.read_csv(f"/Users/oliveirads/Documents/projects/chimerate/ChimeraTE_v2/stringtie/dari_head_TE-initiated_final.tsv", sep ="\t")

for chim_out in x:
    ## Create bed for TE positions
    chim_transcripts = import_chimeras(chim_out)
    chim_transcripts["scaf"] = chim_transcripts["TE_pos"].str.split(":").str[0]
    chim_transcripts["start"] = chim_transcripts["TE_pos"].str.split(":").str[1].str.split("-").str[0]
    chim_transcripts["end"] = chim_transcripts["TE_pos"].str.split(":").str[1].str.split("-").str[1]
    bed_TEs = chim_transcripts[["scaf", "start", "end"]]

    ### Overlap between exons and TEs
    res_overlap = overlap_fraction_intersection(bed_TEs, bed_exons, 0.8)
    
    for row in chimTE_out_df.itertuples():
        scaf = row.TE_pos.split(":")[0].strip()
        start_pos = int(row.TE_pos.split(":")[1].split("-")[0].strip())
        end_pos = int(row.TE_pos.split(":")[1].split("-")[1].strip())

        # Filter overlap_df for matching rows
        matched_rows = res_overlap[
            (res_overlap['Chromosome'].str.strip() == scaf) &
            (res_overlap['Start'] == start_pos) &
            (res_overlap['End'] == end_pos)
        ]

        if not matched_rows.empty:
            # Extract all transcript IDs from 'Name' column
            transcript_ids = matched_rows['Name'].tolist()

            ## Get unique gene ID = remove isoform number
            unique_gene_id = list({tid.rsplit('.', 1)[0] for tid in transcript_ids})
            ## Extract total gene TPM
            matched_genes = stringtie_abund[stringtie_abund['Gene ID'].isin(unique_gene_id)]
            total_gene_TPM = matched_genes['TPM'].iloc[0] if not matched_genes.empty else 0


            # Filter stringtie_TPM_df for matching transcriptIDs
            matched_TPMs = stringtie_TPM_df[
                stringtie_TPM_df['transcriptID'].isin(transcript_ids)]

            ## Compute chimeric isoform TPM
            chimeric_TPM = matched_TPMs['TPM'].sum()

            ## Compute contribution in %
            contribution = (float(chimeric_TPM) * 100 / total_gene_TPM)
            # print(f'{unique_gene_id[0]}\t{total_gene_TPM:.4f}\t{chimeric_TPM:.4f}\t{contribution:.4f}')
            # Store results
            gene_ids.append(unique_gene_id[0])
            total_gene_TPMs.append(f'{total_gene_TPM:.4f}')
            chimeric_TPMs.append(f'{chimeric_TPM:.4f}')
            contributions.append(f'{contribution:.4f}')
        else:
            # Store NaNs if no match
            gene_ids.append(None)
            total_gene_TPMs.append(0)
            chimeric_TPMs.append(0)
            contributions.append(0)


chimTE_out_df['transcript_id'] = gene_ids
chimTE_out_df['total_gene_TPM'] = total_gene_TPMs
chimTE_out_df['chimeric_TPM'] = chimeric_TPMs
chimTE_out_df['contribution_percent'] = contributions

chimTE_out_df.to_csv("chimTE_out_with_TPM.csv", index=False)
# chim_transcripts["start"] = chim_transcripts["TE_pos"].str.split(":").str[1].str.split("-").str[0]
# chim_transcripts["end"] = chim_transcripts["TE_pos"].str.split(":").str[1].str.split("-").str[1]



# TPM_gene = (
#     TPM_stringtie
#     .query('feature == "exon"')
#     [["scaf", "start", "end", "geneID", "strand", "strand"]]
#     .drop_duplicates()
# )

# Dictionary where key = geneID and value = dataframe of that gene
# TPM_gene_dict = {
#     gene_id: all_gene_IDs for gene_id, all_gene_IDs in TPM_gene.groupby("geneID")
# }

# print(TPM_gene_dict)
# for gene_id in all_gene_IDs:
#     TPM_gene = TPM_stringtie.query('geneID == @gene_id and feature == "exon"')[["scaf", "start", "end", "geneID", "strand", "strand"]].drop_duplicates()
    
    
    
#     print(TPM_gene)
#     if TPM_gene.item() < 1:
#         print(f'{gene_id}\t{round(TPM_gene.item(), 4)}')






### Get all geneIDs
# all_gene_IDs = chim_transcripts["gene_id"]#.drop_duplicates().to_list()
# gene_ids = all_gene_IDs["gene_id"]

# Filter TPM_stringtie
# filtered_df = TPM_stringtie[
#     (TPM_stringtie["feature"] == "exon") &
#     (TPM_stringtie["geneID"].isin(all_gene_IDs))]

# bed_exons = filtered_df[["scaf", "start", "end", "geneID", "strand", "strand"]]
# print(bed_exons)