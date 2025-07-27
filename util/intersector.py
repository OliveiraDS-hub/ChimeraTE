import pyranges as pr
import pandas as pd
import subprocess
import numpy as np
import tempfile
import os

def intersection_any_bp(bed1, bed2, out=False):
    ### Load BED files
    if isinstance(bed1, pd.DataFrame):
        if "Start_b" in bed1.columns:
            bed1.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
        bed1_pr = pr.PyRanges(bed1)
    else:
        bed1_pr = pr.read_bed(bed1)

    if isinstance(bed2, pd.DataFrame):
        if "Start_b" in bed2.columns:
            bed2.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
        bed2_pr = pr.PyRanges(bed2)
    else:
        bed2_pr = pr.read_bed(bed2)

    ### Perform intersection: get regions in bed2 that overlap bed1
    overlap = bed1_pr.join(bed2_pr)

    ### Select only BED6 columns BEFORE dropping duplicates
    if overlap.empty:
        return None
    else:
        bed6_df = overlap.df[["Chromosome", "Start", "End", "Name", "Score", "Strand"]]
        bed6_df = bed6_df.drop_duplicates()

        # ID_column = bed6_df.iloc[:, [3]].drop_duplicates().rename(columns={3: "Name"})
        if out is False:
            return bed6_df
        else:
            bed6_df.to_csv(out, header=False, index=False, sep="\t")
            return bed6_df
        
def intersection_any_bp_chunked(bed1, bed2, out_file, chunk_size=1000000):
    if not isinstance(bed1, pd.DataFrame):
        bed1 = pr.read_bed(bed1).df

    if not isinstance(bed2, pd.DataFrame):
        bed2 = pr.read_bed(bed2).df

    result_chunks = []
    for i in range(0, len(bed2), chunk_size):
        chunk = bed2.iloc[i:i+chunk_size]
        chunk_pr = pr.PyRanges(chunk)
        overlap = bed1.join(pr.PyRanges(chunk_pr))
        result_chunks.append(overlap.df)

    pd.concat(result_chunks, ignore_index=True).drop_duplicates().to_csv(out_file, sep ="\t", index=False)

def intersection_any_bp_mbed(bed1, bed2, out=False):
    # Load BED files
    bed1_pr = pr.read_bed(bed1)
    bed2_pr = pr.read_bed(bed2) 

    # Perform intersection: get regions in bed2 that overlap bed1
    overlap = bed1_pr.join(bed2_pr)

    # Drop duplicates based on BED6 fields
    unique = overlap.df.drop_duplicates()
    return unique

def intersection_any_bp_bedtools(bed1, bed2, out_file):
    # Run bedtools and capture output
    with open(out_file, "w") as out_f:
        subprocess.run(
            ["bedtools", "intersect", "-a", bed1, "-b", bed2],
            stdout=out_f, stderr=subprocess.PIPE, check=True)

def intersection_any_bp_bedtools_tmp(bed1, bed2, out_file, tmp_dir):
    subprocess.run(f'bedtools intersect -a {bed1} -b {bed2} > {tmp_dir}/itsect_TEs.bed', shell=True)

    # Step 3: Deduplicate by read ID (column 4)
    # seen_ids = set()
    # with open(f'{tmp_dir}/itsect_TEs.bed', "r") as tmp_in, open(out_file, "w") as out_f:
    #     for line in tmp_in:
    #         cols = line.strip().split("\t")
    #         if len(cols) >= 4:
    #             read_id = cols[3]
    #             if read_id not in seen_ids:
    #                 seen_ids.add(read_id)
    #                 out_f.write(line)

    # # Step 4: Clean up temp file
    # os.remove(tmp_path)

def get_read_ids_from_bedtools2(bed1, bed2):
    cmd = ["bedtools", "intersect", "-a", bed1, "-b", bed2, ">", ]

    # Run the command and capture stdout in memory
    subprocess.run(cmd, shell=True)

    data = np.genfromtxt('your_file.csv', delimiter=',', skip_header=1)
    column_c = data[:, 2]


    # # Process each line from output, extracting column 4
    # ids = []
    # for line in result.stdout.splitlines():
    #     cols = line.split('\t')
    #     ids.append(cols[3])
    # return ids

def overlap_fraction_intersection(bed1_file, bed2_file, fraction=1.0):
    """
    Return DataFrame with 11 columns:
    - 6 columns from bed1 (read)
    - 5 columns from bed2 (range) with suffix '_b'
    Keep only rows where overlap_fraction of bed1 >= `fraction`.
    """
    if isinstance(bed1_file, pd.DataFrame):
        bed1 = bed1_file.copy()
        bed1.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    else:
        try:
            # Try to load with standard BED6 column names
            bed1 = pd.read_csv(bed1_file, sep="\t", header=None,
                            names=["Chromosome", "Start", "End", "Name", "Score", "Strand"])
        except ValueError:
            # Fallback: load with default column names or accept DataFrame directly
            bed1 = pd.read_csv(bed2_file, sep="\t", header=None,
                            names=["Chromosome", "Start_b", "End_b", "Name_b", "Score_b", "Strand_b"])
            bed1.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]

    if isinstance(bed2_file, pd.DataFrame):
        bed2 = bed2_file.copy()
        bed2.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    else:
        try:
            # Try to load with standard BED6 column names
            bed2 = pd.read_csv(bed2_file, sep="\t", header=None,
                            names=["Chromosome", "Start", "End", "Name", "Score", "Strand"])
        except ValueError:
            # Fallback: load with default column names or accept DataFrame directly
            bed2 = pd.read_csv(bed2_file, sep="\t", header=None,
                            names=["Chromosome", "Start_b", "End_b", "Name_b", "Score_b", "Strand_b"])
            bed2.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]

    # Ensure Start/End are integers
    bed1["Start"] = bed1["Start"].astype(int)
    bed1["End"] = bed1["End"].astype(int)
    bed2["Start"] = bed2["Start"].astype(int)
    bed2["End"] = bed2["End"].astype(int)

    # Merge on Chromosome
    merged = bed1.merge(bed2, on="Chromosome", suffixes=("", "_b"))

    # Compute overlap
    start = np.maximum(merged["Start"].values, merged["Start_b"].values)
    end = np.minimum(merged["End"].values, merged["End_b"].values)
    overlap_len = np.clip(end - start, a_min=0, a_max=None)
    merged["overlap_len"] = (overlap_len - start).clip(lower=0)

    # Calculate fraction of bed1 covered by bed2
    merged["read_len"] = merged["End"] - merged["Start"]
    merged["overlap_fraction"] = merged["overlap_len"] / merged["read_len"]

    # Filter rows where overlap_fraction ≥ fraction
    qualified = merged[merged["overlap_fraction"] >= fraction]

    # Select 11 columns: 6 from bed1 + 5 from bed2
    cols = ["Chromosome", "Start", "End", "Name", "Score", "Strand",
            "Start_b", "End_b", "Name_b", "Score_b", "Strand_b"]
    result = qualified[cols].drop_duplicates()

    return result.reset_index(drop=True)


def overlap_fraction_intersection_pyranges(bed1_df, bed2_df, fraction=1.0):
    """
    Intersect bed1 and bed2 using PyRanges.
    Keep rows where bed1 ranges are covered >= `fraction` by bed2 ranges.
    Return 11 columns: 6 from bed1 + 5 from bed2 (_b suffix).
    """

    # Ensure proper types for PyRanges
    if isinstance(bed1_df, pd.DataFrame):
        bed1 = bed1_df.copy()
        bed1.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    else:
        try:
            # Try to load with standard BED6 column names
            bed1 = pd.read_csv(bed1_df, sep="\t", header=None,
                            names=["Chromosome", "Start", "End", "Name", "Score", "Strand"])
        except ValueError:
            # Fallback: load with default column names or accept DataFrame directly
            bed1 = pd.read_csv(bed1_df, sep="\t", header=None,
                            names=["Chromosome", "Start_b", "End_b", "Name_b", "Score_b", "Strand_b"])
            bed1.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]

    if isinstance(bed2_df, pd.DataFrame):
        bed2 = bed2_df.copy()
        bed2.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
    else:
        try:
            # Try to load with standard BED6 column names
            bed2 = pd.read_csv(bed2_df, sep="\t", header=None,
                            names=["Chromosome", "Start", "End", "Name", "Score", "Strand"])
        except ValueError:
            # Fallback: load with default column names or accept DataFrame directly
            bed2 = pd.read_csv(bed2_df, sep="\t", header=None,
                            names=["Chromosome", "Start_b", "End_b", "Name_b", "Score_b", "Strand_b"])
            bed2.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]


    for df in (bed1, bed2):
        df["Chromosome"] = df["Chromosome"].astype(str)
        df["Start"] = pd.to_numeric(df["Start"], downcast="integer").astype("int64")
        df["End"] = pd.to_numeric(df["End"], downcast="integer").astype("int64")
    
    # Convert to PyRanges
    pr_bed1 = pr.PyRanges(bed1)
    pr_bed2 = pr.PyRanges(bed2)

    # Perform interval join (very fast)
    joined = pr_bed1.join(pr_bed2, suffix="_b")

    # Back to Pandas for filtering
    df = joined.df

    if not df.empty:
        # Calculate overlap length
        df["overlap_len"] = (df["End"].clip(upper=df["End_b"]) - df["Start"].clip(lower=df["Start_b"])).clip(lower=0)

        # Fraction of bed1 covered
        df["read_len"] = df["End"] - df["Start"]
        df["overlap_fraction"] = df["overlap_len"] / df["read_len"]

        # Filter by fraction
        qualified = df[df["overlap_fraction"] >= fraction]

        # Select 11 columns and drop duplicates
        cols = ["Chromosome", "Start", "End", "Name", "Score", "Strand",
                "Start_b", "End_b", "Name_b", "Score_b", "Strand_b"]
        result = qualified[cols].drop_duplicates()

        return result.reset_index(drop=True)
    else:
        return None






# def overlap_fraction_intersection(a_bed_file, b_bed_file, fraction=1.0):
#     """
#     Intersect A and B using PyRanges, keeping only A intervals where the fraction
#     of A overlapped by B is >= `fraction`.
#     """
#     import pandas as pd

#     def load_bed(bed_file):
#         if isinstance(bed_file, pd.DataFrame):
#             df = bed_file.copy()
#         else:
#             df = pd.read_csv(bed_file, sep="\t", header=None,
#                              names=["Chromosome", "Start", "End", "Name", "Score", "Strand"])
#         # Fix dtypes BEFORE PyRanges
#         df["Chromosome"] = df["Chromosome"].astype(str)
#         df["Start"] = pd.to_numeric(df["Start"], downcast="integer").astype("int64")
#         df["End"] = pd.to_numeric(df["End"], downcast="integer").astype("int64")
#         return pr.PyRanges(df)

#     a = load_bed(a_bed_file)
#     b = load_bed(b_bed_file)

#     joined = a.join(b, suffix="_b", how="left", report_overlap=True)

#     df = joined.df.copy()
#     df["A_length"] = df["End"] - df["Start"]
#     df["overlap_fraction"] = df["Overlap"] / df["A_length"]

#     df = df[df["overlap_fraction"] >= float(fraction)]

#     if df.empty:
#         return None

#     a_cols = [c for c in df.columns if not c.endswith("_b") and c not in {"Overlap", "A_length", "overlap_fraction"}]
#     b_cols = [c for c in df.columns if c.endswith("_b")]

#     return df[a_cols + b_cols].reset_index(drop=True)




# def intersection_any_bp_bedtools2(bed1, bed2, out_file=None):
#     # Use subprocess to run bedtools
#     cmd = ["bedtools", "intersect", "-a", bed1, "-b", bed2]
    
#     # Stream output directly from bedtools
#     process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

#     ids = []

#     if out_file:
#         with open(out_file, "w") as out_f:
#             for line in process.stdout:
#                 out_f.write(line)
#                 cols = line.strip().split('\t')
#                 if len(cols) >= 4:
#                     ids.append(cols[3])
#     else:
#         for line in process.stdout:
#             cols = line.strip().split('\t')
#             if len(cols) >= 4:
#                 ids.append(cols[3])

#     process.stdout.close()
#     process.wait()

#     return ids
#     # if isinstance(bed1, pd.DataFrame):
    #     bed1_file = "temp_bed1.bed"
    #     bed1[["Chromosome", "Start", "End", "Name", "Score", "Strand"]].to_csv(bed1_file, sep="\t", header=False, index=False)
    # else:
    #     bed1_file = bed1

    # if isinstance(bed2, pd.DataFrame):
    #     bed2_file = "temp_bed2.bed"
    #     bed2[["Chromosome", "Start", "End", "Name", "Score", "Strand"]].to_csv(bed2_file, sep="\t", header=False, index=False)
    # else:
    #     bed2_file = bed2

    # if out is False:
    #     return pd.read_csv(output_file, sep="\t", header=None, names=["Chromosome", "Start", "End", "Name", "Score", "Strand"])


# def overlap_fraction_intersection(a_bed_file, b_bed_file, fraction=1.0):
#     """
#     Intersect A and B using pyranges, keeping only A intervals where the fraction
#     of A overlapped by B is >= `fraction`.

#     Equivalent to: bedtools intersect -a A -b B -wa -wb -f <fraction>
#     """
#     # Read BED files
#     if isinstance(a_bed_file, pd.DataFrame):
#         if "Start_b" in a_bed_file.columns:
#             a_bed_file.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
#         a = pr.PyRanges(a_bed_file)
#     else:
#         a = pr.read_bed(a_bed_file)

#     if isinstance(b_bed_file, pd.DataFrame):
#         if "Start_b" in b_bed_file.columns:
#             b_bed_file.columns = ["Chromosome", "Start", "End", "Name", "Score", "Strand"]
#         b = pr.PyRanges(b_bed_file)
#     else:
#         b = pr.read_bed(b_bed_file)

#     # Ensure int type for start and end
#     a.df["Start"] = a.df["Start"].astype("int64")
#     a.df["End"] = a.df["End"].astype("int64")
#     b.df["Start"] = b.df["Start"].astype("int64")
#     b.df["End"] = b.df["End"].astype("int64")

#     # Join and compute overlap
#     joined = a.join(b, suffix="_b", how="left", report_overlap=True)

#     df = joined.df.copy()
#     df["A_length"] = df["End"] - df["Start"]
#     df["overlap_fraction"] = df["Overlap"] / df["A_length"]

#     # Filter for required fraction
#     df = df[df["overlap_fraction"] >= float(fraction)]

#     if df.empty:
#         return None
    
#     #  Include all original A and B columns
#     a_cols = [c for c in df.columns if not c.endswith("_b") and c not in {"Overlap", "A_length", "overlap_fraction"}]
#     b_cols = [c for c in df.columns if c.endswith("_b")]

#     return df[a_cols + b_cols].reset_index(drop=True)