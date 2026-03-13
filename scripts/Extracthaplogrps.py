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
