import pandas as pd
import numpy as np

#######################################
# 1. LOAD DATA
#######################################

# Load the AADR annotations file, selecting only relevant columns
columns = [
    'Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]',
    'Group ID', 'Political Entity', 'Lat.', 'Long.',
    'Molecular Sex', 'Y haplogroup (manual curation in terminal mutation format)', 
    'Y haplogroup (manual curation in ISOGG format)', 'mtDNA haplogroup if >2x or published'
     ]

df = pd.read_csv(
    "/home/inf-41-2025/BINP29/Popgenetics/Project/AADR_Annotations_2025.tsv",
    sep="\t",
    usecols=columns,
    engine="python",
    on_bad_lines="skip")

print(df.head(2))

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

def create_frequency_table(df, hap_column, outfile, include_sex=False):
    """Create frequency tables and basal-level haplists."""
    df = df.copy()
    df["basal"] = df[hap_column].str.extract(r'^([A-Z])')

    # Frequency table
    if include_sex:
        freq = df.groupby(["Ancient pop","Sex","basal"]).size().unstack(fill_value=0).reset_index()
        freq["total"] = freq.select_dtypes(include="number").sum(axis=1)
        meta = df.groupby(["Ancient pop","Sex"]).agg({
            "Age": "mean",
            "Country": "first",
            "Lat": "mean",
            "Long": "mean"
        }).reset_index()
        meta["Age"] = meta["Age"].round().astype(int)
        final = meta.merge(freq, on=["Ancient pop","Sex"], how="inner")
        final = final[~final["Sex"].str.startswith(("c", "U"), na=False)]
    else:
        freq = df.groupby(["Ancient pop","basal"]).size().unstack(fill_value=0).reset_index()
        freq["total"] = freq.select_dtypes(include="number").sum(axis=1)
        meta = df.groupby("Ancient pop").agg({
            "Age": "mean",
            "Country": "first",
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

df_yter = clean_haplogroup(df, "Y_Haplogroup")
len(df_yter)
print(df_yter.head(2))
df_yter_freq = create_frequency_table(df_yter, "Y_Haplogroup",
                       "/home/inf-41-2025/BINP29/Popgenetics/y_hap_ter_freq.tsv")
len(df_yter_freq)
print(df_yter_freq.head(2))
