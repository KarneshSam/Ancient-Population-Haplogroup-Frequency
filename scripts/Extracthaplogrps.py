import pandas as pd
import numpy as np
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

# Rename columns to more user-friendly names
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

print(df.head(2))

#######################################
# 3. CLEAN COORDINATES AND ANCIENT POPULATION SAMPLES
#######################################

# Remove rows with missing coordinates in Long and Lat
df = df[df["Lat"].replace("..", np.nan).notna() &
        df["Long"].replace("..", np.nan).notna()]
# Convert Lat and Long to numeric, coercing errors to NaN
df["Lat"] = pd.to_numeric(df["Lat"], errors="coerce")
df["Long"] = pd.to_numeric(df["Long"], errors="coerce")
# Drop rows where Lat or Long is still NaN after conversion
df = df.dropna(subset=["Lat", "Long"])

print(df.head(2))
print(df[["Lat", "Long"]].head(2))

# Convert Lat and Long to float, replacing commas with dots
df["Lat"] = df["Lat"].str.replace(",", ".").astype(float)
df["Long"] = df["Long"].str.replace(",", ".").astype(float)

print(df[["Lat", "Long"]].head(2))

# Keep rows which are Ancient pop
df = df[df["Age"] != 0]
len(df)

#######################################
# 4. FILTER N/A VALUES IN HAPLOGROUP COLUMNS
#######################################

# Define a function to clean haplogroup columns
# This function removes rows where the specified column has missing values (represented as "..") 
# or values that start with "n/a", "na", "NO", or "not published"
def clean_haplogroup(df, column):
    """Remove missing or invalid haplogroups"""
    temp = df[df[column].replace("..", np.nan).notna()]
    temp = temp[~temp[column].str.startswith(("n/a", "na", "NO", "not published"), na=False)]
    return temp

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
        haplists_file = outfile.replace(".tsv", "_haplists_by_sex.tsv")
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
    return final

#######################################
# 6. Y HAPLOGROUP TERMINAL
#######################################
df_yter = clean_haplogroup(df, "Y_Haplogroup")
len(df_yter)
print(df_yter.head(2))
df_yter_freq = create_frequency_table(df_yter, "Y_Haplogroup",
                       "/home/inf-41-2025/BINP29/Popgenetics/y_hap_ter_freq.tsv")
len(df_yter_freq)
print(df_yter_freq.head(2))

#######################################
# 7. Y HAPLOGROUP ISOGG
#######################################
df_yisogg = clean_haplogroup(df, "Y_Haplogroup_ISOGG")
len(df_yisogg)
print(df_yisogg.head(2))
df_yisogg_freq = create_frequency_table(df_yisogg, "Y_Haplogroup_ISOGG",
                       "/home/inf-41-2025/BINP29/Popgenetics/y_hap_isogg_freq.tsv")
len(df_yisogg_freq)
print(df_yisogg_freq.head(2))   

#######################################
# 8. mtDNA HAPLOGROUP  
#######################################
df_mt = clean_haplogroup(df, "mtDNA_Haplogroup")
len(df_mt)
print(df_mt.head(2))
df_mt_freq = create_frequency_table(df_mt, "mtDNA_Haplogroup",
                       "/home/inf-41-2025/BINP29/Popgenetics/mt_hap_freq.tsv", include_sex=True)
len(df_mt_freq)
print(df_mt_freq.head(2))