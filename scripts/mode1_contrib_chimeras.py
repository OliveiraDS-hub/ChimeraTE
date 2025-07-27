import pyranges as pr
import pandas as pd
import numpy as np
import subprocess
import warnings
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
import argparse
import shutil
import os
import re

parser = argparse.ArgumentParser(description='')
parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')
required.add_argument('--input', help='Paired-end files and their respective group/replicate', required=True, type=str, metavar = "")
required.add_argument('--index', help='Directory STAR index', required=True, type=str, metavar = "")
required.add_argument('--output', help='prefix for output files', required=True, type=str, metavar = "")
required.add_argument('--gff', help='gff with gene annotation', required=True, type=str, metavar = "")
required.add_argument('--dir_chimerate', help='directory with chimeraTE output tsv files', required=True, type=str, metavar= "")
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
    stringtie_transcripts = stringtie_df[stringtie_df["feature"] == "transcript"].copy()

    ##stringtie_transcripts['TPM'] = stringtie_transcripts['ID'].str.extract(r'TPM "([^"]+)"').astype(float)
    stringtie_transcripts.loc[:, 'TPM'] = (
        stringtie_transcripts['ID'].str.extract(r'TPM "([^"]+)"').astype(float))
    
    stringtie_transcripts['ref_gene'] = stringtie_transcripts['ID'].str.extract(r'ref_gene_name "([^"]+)"').astype(str)

    # Convert string "nan" to real NaN
    # stringtie_transcripts['ref_gene'].replace("nan", pd.NA, inplace=True)
    stringtie_transcripts['ref_gene'] = stringtie_transcripts['ref_gene'].replace("nan", pd.NA)

    # Extract gene_id from transcriptID (STRG.22.1 → STRG.22)
    stringtie_transcripts['gene_id'] = stringtie_transcripts['transcriptID'].str.extract(r'(STRG\.\d+)')

    # Fill NaNs in ref_gene per gene_id
    stringtie_transcripts['ref_gene'] = (
        stringtie_transcripts.groupby('gene_id')['ref_gene']
        .transform(lambda x: x.ffill().bfill())
    )

    # Keep only ref_gene, transcriptID, TPM
    tpm_transcripts = stringtie_transcripts[['ref_gene', 'transcriptID', 'TPM']]
    tpm_transcripts.to_csv(f"{out_prefix}_stringtie_TPM_values.csv", index=False)

    ## extract bed exons
    stringtie_exons = stringtie_df[stringtie_df["feature"] == "exon"]
    stringtie_exons.loc[:, 'TPM'] = (
        stringtie_exons['ID'].str.extract(r'TPM "([^"]+)"').astype(float))
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
            # a_bed_file.to_csv("TEs.bed", sep = "\t", index = False, header=None)
        a = pr.PyRanges(a_bed_file)
    else:
        a = pr.read_bed(a_bed_file)
    
    if isinstance(b_bed_file, pd.DataFrame):
        if len(b_bed_file.columns) == 6:
            b_bed_file.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
            # b_bed_file.to_csv("exons.bed", sep = "\t", index = False, header=None)
        b = pr.PyRanges(b_bed_file)
    else:
        b = pr.read_bed(b_bed_file)

    # Join and compute overlap
    joined = a.join(b, suffix="_b", how="left", report_overlap=True)

    df = joined.df.copy()
    df["TE_length"] = df["End"] - df["Start"]
    df["overlap_fraction"] = df["Overlap"] / df["TE_length"]
    
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
    print(f"\nMerged files have been found!\tSkipping...")

### STAR alignment
if not os.path.exists(f"{args.output}_Aligned.sortedByCoord.out.bam"):
    subprocess.run(f'STAR --genomeDir {args.index} --runThreadN 16 --readFilesCommand zcat --readFilesIn merged_R1.fastq.gz merged_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --outFileNamePrefix {args.output}_ --outTmpDir STARtmp', shell=True)
else:
    print("\nAlignment bam has been found!")


### Stringtie expression
if not os.path.exists(f'{args.output}_gene_sample.gtf'):
    subprocess.run(f'stringtie {args.output}_Aligned.sortedByCoord.out.bam -G {args.gff} -o {args.output}_gene_sample.gtf -A {args.output}_gene_abundances.tsv -B -p 8', shell=True)
else:
    print("\nstringtie output has been found!\tskipping...")

## Create bed for exon positions
bed_exons = getTPM(f"{args.output}_gene_sample.gtf", str(args.output))

## Import stringtie outputs
stringtie_TPM_df = pd.read_csv(f"{args.output}_stringtie_TPM_values.csv")
stringtie_abund = pd.read_csv(f"{args.output}_gene_abundances.tsv", usecols=[0,8], sep ="\t")

chimTE_outputs = ['TE-initiated_final.tsv',
                  'TE-exonized_embedded_final.tsv',
                  'TE-exonized_intronic_final.tsv',
                  'TE-exonized_overlapped_final.tsv',
                  'TE-terminated_final.tsv']

for chim_tsv in chimTE_outputs:
    chim_type = re.sub("_final.tsv", "", chim_tsv)
    print(f'Analyzing contribution of chimeric transcripts from {chim_type}...', end = "\t", flush=True)

    ### Store the results
    gene_ids = []
    total_gene_TPMs = []
    chimeric_TPMs = []
    contributions = []

    ## Import dataframe
    chim_out = pd.read_csv(f'{args.dir_chimerate}/{chim_tsv}', sep ="\t")

    ## Create bed for TE positions
    chim_transcripts = import_chimeras(f'{args.dir_chimerate}/{chim_tsv}')
    chim_transcripts["scaf"] = chim_transcripts["TE_pos"].str.split(":").str[0]
    chim_transcripts["start"] = chim_transcripts["TE_pos"].str.split(":").str[1].str.split("-").str[0]
    chim_transcripts["end"] = chim_transcripts["TE_pos"].str.split(":").str[1].str.split("-").str[1]
    bed_TEs = chim_transcripts[["scaf", "start", "end"]]

    ### Overlap between exons and TEs
    res_overlap = overlap_fraction_intersection(bed_TEs, bed_exons, 0.8)

    for row in chim_out.itertuples():
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
            unique_gene_ids = list({tid.rsplit('.', 1)[0] for tid in transcript_ids})
            
            for uniq_gene in unique_gene_ids:
                if uniq_gene in gene_ids:
                    continue
                ## Extract total gene TPM
                matched_genes = stringtie_abund[stringtie_abund['Gene ID'].isin([uniq_gene])]
                total_gene_TPM = matched_genes['TPM'].iloc[0] if not matched_genes.empty else 0

                # subset transcript ids only for the this specific uniq gene ID
                specific_transcripts = [item for item in transcript_ids if item.startswith(uniq_gene)]

                # Filter stringtie_TPM_df for matching transcriptIDs
                matched_TPMs = stringtie_TPM_df[
                    stringtie_TPM_df['transcriptID'].isin(specific_transcripts)]

                ## Compute chimeric isoform TPM
                chimeric_TPM = matched_TPMs['TPM'].sum()

                ## Compute contribution in %
                contribution = (float(chimeric_TPM) * 100 / total_gene_TPM)

                # if "STRG.11661" in uniq_gene:
                #     # print(specific_transcripts)
                #     # print(matched_genes)
                #     print(f'{uniq_gene}\t{total_gene_TPM}')#\t{chimeric_TPM}\t{contribution}')
                # elif "STRG.11662" in uniq_gene:
                #     # print(specific_transcripts)
                #     print(f'{uniq_gene}\t{total_gene_TPM}')

            ### Store results
            gene_ids.append(uniq_gene)
            total_gene_TPMs.append(f'{total_gene_TPM:.4f}')
            chimeric_TPMs.append(f'{chimeric_TPM:.4f}')
            contributions.append(f'{contribution:.4f}')
        else:
            # Store NaNs if no match
            gene_ids.append(None)
            total_gene_TPMs.append(0)
            chimeric_TPMs.append(0)
            contributions.append(0)

    # ### Export results as tsv
    chim_out['transcript_id'] = gene_ids
    chim_out['total_gene_TPM'] = total_gene_TPMs
    chim_out['chimeric_TPM'] = chimeric_TPMs
    chim_out['contribution_percent'] = contributions

    chim_out.to_csv(f"{args.output}_TPM_{chim_tsv}", sep="\t", index=False)
    print("Done!")