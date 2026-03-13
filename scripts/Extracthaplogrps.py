#!/usr/bin/env python3

"""
Extracthaplogrps.py (script's name)

Description:
      This program processes the AADR annotations dataset and extracts basal 
      haplogroup frequencies for ancient populations. The script reads the 
      annotations file, filters relevant columns, cleans coordinate and haplogroup 
      information, and produces frequency tables for Y-chromosome and mitochondrial 
      haplogroups.
      
      The main goal is to summarize haplogroup distributions across ancient
      populations. The program generates frequency tables for:
            1. Y haplogroups in terminal mutation format
            2. Y haplogroups in ISOGG format
            3. mtDNA haplogroups (calculated separately by sex)

      Additionally, the script produces haplogroup list files (haplists) that
      store the individual haplogroups belonging to each basal haplogroup group
      for each ancient population. These lists allow inspection of the detailed
      haplogroup composition within the summarized frequency tables.

      The script is designed to handle missing values, malformed entries, and
      coordinate inconsistencies commonly found in large annotation datasets.
      It also ensures robustness through error handling for missing files,
      invalid inputs, and file overwriting.

      If the user specifies output filenames, those will be used. Otherwise,
      the program automatically generates default filenames inside the output
      directory. The program also supports forced overwriting of existing files.

List of Functions:
      1. load_data – Reads the AADR annotations file and extracts only the
                     relevant columns needed for analysis.

      2. rename_columns – Renames the original dataset column names to simpler,
                          standardized names used throughout the script.

      3. clean_coordinates – Cleans latitude and longitude values, removes rows
                             with missing coordinates, converts coordinates to
                             numeric values, and filters for ancient samples.

      4. clean_haplogroup – Filters invalid or missing haplogroup values from
                            the specified haplogroup column.

      5. create_frequency_table – Generates haplogroup frequency tables based
                                  on basal haplogroup categories and produces
                                  haplogroup list (haplist) files.

      6. main – Controls the overall workflow including argument parsing,
                 file validation, data processing, and result generation.

List of non-standard modules:
      1. pandas – Used for data manipulation and dataframe operations.
      2. numpy – Used for handling missing values and numeric operations.

Procedure:
      1. The program reads the AADR annotations file using pandas.
      2. The selected columns are renamed to simplified column names.
      3. Geographic coordinates are cleaned and rows with missing values
         are removed and kept only ancient population samples.
      4. Haplogroup columns are filtered to remove invalid or missing entries.
      5. Separate dataframes are created for:
            • Y haplogroup (terminal mutation format)
            • Y haplogroup (ISOGG format)
            • mtDNA haplogroup
      6. Basal haplogroups (single-letter classification) are extracted.
      7. Frequency tables are computed for each ancient population.
      8. Haplogroup lists are generated to retain detailed haplogroup
         information.
      9. Results are written as tab-separated (.tsv) files.

Input:
      1. AADR annotations file (TSV format) containing metadata for ancient
         individuals including haplogroup assignments and geographic data.

Output:
      1. Y haplogroup terminal frequency table (.tsv)
      2. Y haplogroup ISOGG frequency table (.tsv)
      3. mtDNA haplogroup frequency table by sex (.tsv)
      4. Haplogroup list files corresponding to each frequency table.

Usage: python Extracthaplogrps.py -i input_file -o output_dir

Get Help: python Extracthaplogrps.py -h

    Options:
       --yter OUTPUT_FILE: Specify output file for Y haplogroup terminal frequencies.
       --yisogg OUTPUT_FILE: Specify output file for Y haplogroup ISOGG frequencies.
       --mt OUTPUT_FILE: Specify output file for mtDNA haplogroup frequencies.
       -f, --force: Overwrite output files if they already exist.

Version: 1.00
Date: 2026-03-13
Name: Karnesh Sampath

"""

# Pandas and NumPy for data manipulation
# argparse for command-line argument parsing 
# os and sys for file handling and system operations
import pandas as pd
import numpy as np
import argparse
import os
import sys

#######################################
# 1. LOAD DATA
#######################################

# Load the AADR annotations file, selecting only relevant columns
def load_data(input_file):
    """Load the AADR annotations file with selected columns."""
    columns = [
        'Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]',
        'Group ID', 'Political Entity', 'Lat.', 'Long.',
        'Molecular Sex', 'Y haplogroup (manual curation in terminal mutation format)', 
        'Y haplogroup (manual curation in ISOGG format)', 'mtDNA haplogroup if >2x or published'
    ]
    # Use try-except to catch potential errors during file loading
    try:
        df = pd.read_csv(
            input_file,
            sep="\t",
            usecols=columns,
            engine="python",
            on_bad_lines="skip")
    except Exception as e: 
        sys.exit(f"Error loading data: {e}")
    return df
#print(df.head(2))

#######################################
# 2. RENAME COLUMNS
#######################################

# Rename columns to more concise names for easier handling
def rename_columns(df):
    """Rename columns to more concise names for easier handling."""
    df.rename(columns={
      'Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]': 'Age',
      'Group ID': 'Ancient pop',
      'Political Entity': 'Country',
      'Lat.': 'Lat',  
      'Long.': 'Long',
      'Molecular Sex': 'Sex',
      'Y haplogroup (manual curation in terminal mutation format)': 'Y_Haplogroup',
      'Y haplogroup (manual curation in ISOGG format)': 'Y_Haplogroup_ISOGG',
      'mtDNA haplogroup if >2x or published': 'mtDNA_Haplogroup'
      }, inplace=True)
    
    return df

#######################################
# 3. CLEAN COORDINATES AND ANCIENT POPULATION SAMPLES
#######################################

# Define a function to clean the Lat and Long columns by removing entire rows with missing values
def clean_coordinates(df):
    """Clean the Lat and Long columns by removing rows with missing values and converting to numeric."""
    df = df.copy()
    # Remove rows with missing coordinates in Long and Lat
    df = df[df["Lat"].replace("..", np.nan).notna() &
        df["Long"].replace("..", np.nan).notna()]
    
    # Convert Lat and Long to float, replacing commas with dots
    df["Lat"] = df["Lat"].astype(str).str.replace(",", ".")
    df["Long"] = df["Long"].astype(str).str.replace(",", ".")

    # Convert Lat and Long to numeric, coercing errors to NaN
    df["Lat"] = pd.to_numeric(df["Lat"], errors="coerce")
    df["Long"] = pd.to_numeric(df["Long"], errors="coerce")
   
    # Drop rows where Lat or Long is still NaN after conversion
    df = df.dropna(subset=["Lat", "Long"])

    # Keep rows which are Ancient pop
    df = df[df["Age"] != 0]

    return df

#######################################
# 4. FILTER N/A VALUES IN HAPLOGROUP COLUMNS
#######################################

# Define a function to clean haplogroup columns
# This function removes rows where the specified column has missing values (represented as "..") 
# or values that start with "n/a", "na", "NO", or "not published"
def clean_haplogroup(df, column):
    """Remove missing or invalid haplogroups"""
    df = df.copy()
    df = df[df[column].replace("..", np.nan).notna()]
    df = df[~df[column].str.startswith(("n/a", "na", "NO", "not published"), na=False)]
    return df

#######################################
# 5. CREATE SEPARATE DATAFRAMES FOR EACH HAPLOGROUP TYPE
#######################################

# Create separate dataframes for each haplogroup type
# This creates a basal (single letter) haplogroup column for each type, which will be used for frequency tables
def create_frequency_table(df, hap_column, outfile, include_sex=False):
    """Create frequency tables and basal-level haplists."""
    df = df.copy()
    df["basal"] = df[hap_column].str.extract(r'^([A-Z])')

    # Basal haplists
    # for sex-specific haplogroup list
    if include_sex:
        # Extract the haplogroups for each Ancient pop with repect to the basal haplogroup
        hap_lists = (
            df.groupby(["Ancient pop", "Sex", "basal"])[hap_column]
              .apply(lambda x: ",".join(x))
              .unstack()
              ).reset_index()
        hap_lists = hap_lists[~hap_lists["Sex"].str.startswith(("c", "U"), na=False)] # Remove rows where Sex starts with "c" or "U" (child or unknown)
        haplists_file = outfile.replace(".tsv", "_haplists.tsv")
    # for non sex-specific haplogroup list    
    else:
        hap_lists = (
            df.groupby(["Ancient pop", "basal"])[hap_column]
              .apply(lambda x: ",".join(x))
              .unstack()
              ).reset_index()
        haplists_file = outfile.replace(".tsv", "_haplists.tsv")
    # Save the haplogroup lists to a separate file
    hap_lists.to_csv(haplists_file, sep="\t", index=False)

    # Frequency table
    if include_sex:
        freq = df.groupby(["Ancient pop","Sex","basal"]).size().unstack(fill_value=0).reset_index()
        freq["total"] = freq.select_dtypes(include="number").sum(axis=1)      # total count for each Ancient pop; includes only the numeric columns (the haplogroups)
        # Metadata table with mean age, first country, mean lat and long for each Ancient pop
        meta = df.groupby(["Ancient pop","Sex"]).agg({
            "Age": "mean",
            "Country": "first",
            "Lat": "mean",
            "Long": "mean"
            }).reset_index()
        meta["Age"] = meta["Age"].round().astype(int)
        final = meta.merge(freq, on=["Ancient pop","Sex"], how="inner")
        final = final[~final["Sex"].str.startswith(("c", "U"), na=False)] # Remove rows where Sex starts with "c" or "U" (child or unknown)
    else:
        freq = df.groupby(["Ancient pop","basal"]).size().unstack(fill_value=0).reset_index()
        freq["total"] = freq.select_dtypes(include="number").sum(axis=1)
        # Metadata table with mean age, first country, mean lat and long for each Ancient pop
        meta = df.groupby("Ancient pop").agg({
            "Age": "mean",
            "Country": "first",
            "Sex": "first",
            "Lat": "mean",
            "Long": "mean"
            }).reset_index()
        meta["Age"] = meta["Age"].round().astype(int)
        final = meta.merge(freq, on="Ancient pop", how="inner")
    
    # Clean Ancient pop names ("_" to " ")
    if "Ancient pop" in final.columns:
        final["Ancient pop"] = final["Ancient pop"].astype(str).str.replace("_", " ")

    final.to_csv(outfile, sep="\t", index=False)
    print(f"Saved frequency table: {outfile}")
    print(f"Saved haplist table: {haplists_file}")
    return final

#######################################
# MAIN EXECUTION
#######################################

# Define the main function to orchestrate the workflow
def main():
    parser = argparse.ArgumentParser(description="Extract haplogroup frequencies from AADR annotations file")
    # Define command-line arguments for input file, output directory, and optional output files for each haplogroup type
    parser.add_argument(
        "-i", "--input", required=True,
        help="Path to the AADR annotations file (tsv format)")
    parser.add_argument(
        "-o", "--outdir", required=True,
        help="Directory to save the output frequency tables and haplogroup lists")
    parser.add_argument(
        "--yter", 
        help="Output TSV file of Y haplogroup terminal frequencies")
    parser.add_argument(
        "--yisogg", 
        help="Output TSV file of Y haplogroup ISOGG frequencies")
    parser.add_argument(
        "--mt", 
        help="Output TSV file of mtDNA haplogroup frequencies")
    # Add the --force argument
    parser.add_argument(
        "-f", "--force", action="store_true",
        help="Overwrite output files if they already exist")
    # Parse the command-line arguments
    args = parser.parse_args()
    # Mention the input file and output directory
    input_file = args.input
    outdir = args.outdir
    # Set default output file paths if not provided by the user
    yter = args.yter if args.yter else os.path.join(outdir, "y_hap_ter_freq.tsv")
    yisogg = args.yisogg if args.yisogg else os.path.join(outdir, "y_hap_isogg_freq.tsv")
    mt = args.mt if args.mt else os.path.join(outdir, "mt_hap_freq.tsv")
    
    # Check if the input file and output directory exist
    if not os.path.exists(input_file):
        sys.exit(f"Error: Input file '{input_file}' does not exist.")
    if not os.path.exists(outdir):  
        sys.exit(f"Error: Output directory '{outdir}' does not exist.")
    
    # Check if files exist and force overwrite is not set
    for f in [yter, yisogg, mt]:
       if os.path.exists(f) and not args.force:
        sys.exit(f"Error: Output file '{f}' already exists. Use --force to overwrite.")

    # ------- WORKFLOW -------
    # 1. Load data
    df = load_data(input_file)
    # 2. Rename columns
    df = rename_columns(df)
    # 3. Clean coordinates and filter for Ancient pop samples
    df = clean_coordinates(df)
    # 4. Create separate dataframes for each haplogroup type and generate frequency tables
    # Y haplogroup terminal
    df_yter = clean_haplogroup(df, "Y_Haplogroup")
    create_frequency_table(df_yter, "Y_Haplogroup", yter)
    # Y haplogroup ISOGG
    df_yisogg = clean_haplogroup(df, "Y_Haplogroup_ISOGG")
    create_frequency_table(df_yisogg, "Y_Haplogroup_ISOGG", yisogg)
    # mtDNA haplogroup
    df_mt = clean_haplogroup(df, "mtDNA_Haplogroup")
    create_frequency_table(df_mt, "mtDNA_Haplogroup", mt, include_sex=True)
    
##################################
# Run the main function when the script is executed
##################################

if __name__ == "__main__":
    main()
