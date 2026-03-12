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
