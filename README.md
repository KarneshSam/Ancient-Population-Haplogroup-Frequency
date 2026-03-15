# Ancient Population Haplogroup Explorer Pipeline

## Overview

This project implements a two-script workflow for extracting, summarizing, and interactively visualizing ancient population haplogroup distributions from the AADR (Allen Ancient DNA Resource) annotations dataset.

**Script 1 — `Extracthaplogrps.py`** processes the raw AADR annotations file and produces haplogroup frequency tables and subhaplogroup lists for Y-chromosome (terminal and ISOGG formats) and mitochondrial DNA haplogroups across ancient populations.

**Script 2 — `Haplogrpviewer.py`** is an interactive Streamlit web application that loads the outputs from Script 1 and allows users to explore haplogroup distributions geographically and inspect subhaplogroup hierarchies through maps, pie charts, and sunburst diagrams.

---

## Purpose of the Analysis

The main goals of this pipeline are:

* Extract and clean haplogroup assignments from the AADR annotations dataset.
* Summarize Y-chromosome and mitochondrial haplogroup frequencies per ancient population.
* Produce structured output files suitable for downstream visualization and analysis.
* Allow interactive geographic and compositional exploration of haplogroup data through a web interface.

This type of analysis is commonly used in:

* Archaeogenetics and ancient DNA research
* Population genetics and human migration studies
* Comparative genomics across prehistoric populations

---

## Data Description

### Input Data

* **AADR annotations file** (TSV format): Contains metadata for ancient individuals including haplogroup assignments (Y and mtDNA), geographic coordinates, molecular sex, political entity, group ID, and radiocarbon/contextual age estimates.

> *Note:* The input file path and output directory are specified at runtime via command-line arguments. See [Usage](#usage) below.

---

## Workflow

The pipeline is organized into two sequential stages:

### Stage 1 — Haplogroup Extraction (`Extracthaplogrps.py`)

**1. Load Data**
* Purpose: Read the AADR annotations TSV file and select only the relevant columns.
* Input: AADR annotations file (TSV).
* Output: In-memory dataframe with selected columns.
* Tool: `pandas.read_csv` with `usecols`.

**2. Rename Columns**
* Purpose: Standardize verbose column names to concise identifiers used throughout the script.
* Details: Renames columns such as `Y haplogroup (manual curation in terminal mutation format)` to `Y_Haplogroup`.

**3. Clean Coordinates and Filter Ancient Samples**
* Purpose: Remove rows with missing or malformed geographic coordinates and retain only ancient population samples (Age ≠ 0).
* Details: Handles `".."` placeholders, comma-decimal formats, and non-numeric values. Coordinates are converted to float after cleaning.

**4. Clean Haplogroup Columns**
* Purpose: Remove rows with missing or invalid haplogroup entries.
* Details: Filters out values that are `".."`, `"n/a"`, `"na"`, `"NO"`, or `"not published"`.

**5. Generate Frequency Tables and Haplogroup Lists**
* Purpose: Compute basal haplogroup frequencies per ancient population and sex, and export subhaplogroup member lists.
* Output:
  * Frequency TSV: counts per basal haplogroup per population, with metadata (mean age, country, mean lat/long, total count).
  * Haplist TSV: comma-separated lists of all individual haplogroups belonging to each basal haplogroup per population.

Three separate runs are performed for:
* Y haplogroup (terminal mutation format)
* Y haplogroup (ISOGG format)
* mtDNA haplogroup

---

### Stage 2 — Interactive Viewer (`Haplogrpviewer.py`)

**1. Argument Parsing and File Validation**
* Purpose: Accept paths to all six input files (three frequency tables + three haplist files) and validate that they exist before loading.

**2. Load and Cache Datasets**
* Purpose: Load all six TSV files into a nested dictionary of pandas dataframes.
* Tool: `@st.cache_data` to prevent reloading on every user interaction.

**3. Sidebar Filters**
* Purpose: Allow users to interactively filter the dataset by age range, country, and sex.
* Details: Filters are applied using `@st.cache_data`-wrapped function to avoid redundant computation.

**4. Population Table**
* Purpose: Display the filtered dataset as a selectable table, allowing users to click a row to highlight a population.

**5. Geographic Map**
* Purpose: Render an interactive map showing all filtered ancient population locations as clustered markers.
* Tool: `folium` with `MarkerCluster`. 
* Details: Map is cached and only rebuilt when filters change. Live selection state (red marker) is added outside the cache.

**6. Pie Chart**
* Purpose: Display basal haplogroup composition for the selected population.
* Tool: `plotly.express.pie`.

**7. Subhaplogroup Table and Sunburst Diagram**
* Purpose: Display a table of all subhaplogroups per basal haplogroup and a hierarchical sunburst diagram illustrating haplogroup substructure.
* Tool: `plotly.express.sunburst`.

---

## Pipeline Directory Structure

```text
.
├── data/
│   └── AADR_annotations.tsv          # Input AADR annotations file
├── results/
│   ├── y_hap_ter_freq.tsv            # Y terminal frequency table
│   ├── y_hap_ter_freq_haplists.tsv   # Y terminal haplogroup lists
│   ├── y_hap_isogg_freq.tsv          # Y ISOGG frequency table
│   ├── y_hap_isogg_freq_haplists.tsv # Y ISOGG haplogroup lists
│   ├── mt_hap_freq.tsv               # mtDNA frequency table
│   └── mt_hap_freq_haplists.tsv      # mtDNA haplogroup lists
├── scripts/
|   ├── Extracthaplogrps.py           # Script 1: Extraction and frequency tables
|   └── Haplogrpviewer.py             # Script 2: Interactive Streamlit viewer
└── README.md                         # This file
```

---

## Required Software

This workflow depends on the following Python packages:

| Tool | Purpose |
|------|---------|
| `pandas` | Data loading, filtering, and manipulation |
| `numpy` | Handling missing values and numeric operations |
| `argparse` | Command-line argument parsing |
| `streamlit` | Interactive web application framework |
| `folium` | Geographic map visualization |
| `folium.plugins.MarkerCluster` | Marker clustering for map performance |
| `plotly.express` | Pie charts and sunburst diagrams |
| `streamlit_folium` | Embedding Folium maps in Streamlit |

Install all dependencies:

```bash
pip install pandas numpy streamlit folium streamlit-folium plotly
```

---

## Usage

### Script 1 — Extract Haplogroup Frequencies

```bash
python scripts/Extracthaplogrps.py -i input_file -o output_dir
```

**Required Arguments:**

| Argument | Description |
|----------|-------------|
| `-i`, `--input` | Path to the AADR annotations file (TSV format) |
| `-o`, `--outdir` | Directory to save all output files |

**Optional Arguments:**

| Argument | Description |
|----------|-------------|
| `--yter OUTPUT_FILE` | Custom output path for Y terminal frequency table |
| `--yisogg OUTPUT_FILE` | Custom output path for Y ISOGG frequency table |
| `--mt OUTPUT_FILE` | Custom output path for mtDNA frequency table |
| `-f`, `--force` | Overwrite output files if they already exist |

**Get Help:**

```bash
python scripts/Extracthaplogrps.py -h
```

**Example:**

```bash
python scripts/Extracthaplogrps.py -i data/AADR_annotations.tsv -o results/ --force
```

---

### Script 2 — Launch Interactive Viewer

```bash
streamlit run scripts/Haplogrpviewer.py -- \
    --y_term results/y_hap_ter_freq.tsv \
    --y_iso results/y_hap_isogg_freq.tsv \
    --mt results/mt_hap_freq.tsv \
    --y_term_sub results/y_hap_ter_freq_haplists.tsv \
    --y_iso_sub results/y_hap_isogg_freq_haplists.tsv \
    --mt_sub results/mt_hap_freq_haplists.tsv
```

**Required Arguments:**

| Argument | Description |
|----------|-------------|
| `--y_term` | Y haplogroup terminal frequency table (.tsv) |
| `--y_iso` | Y haplogroup ISOGG frequency table (.tsv) |
| `--mt` | mtDNA haplogroup frequency table (.tsv) |
| `--y_term_sub` | Y terminal subhaplogroup list (.tsv) |
| `--y_iso_sub` | Y ISOGG subhaplogroup list (.tsv) |
| `--mt_sub` | mtDNA subhaplogroup list (.tsv) |

**Get Help:**

```bash
streamlit run scripts/Haplogrpviewer.py -- --help
```

---

## Output Files

### From `Extracthaplogrps.py`

| File | Description |
|------|-------------|
| `y_hap_ter_freq.tsv` | Y haplogroup (terminal) frequency table |
| `y_hap_ter_freq_haplists.tsv` | Individual haplogroup members per basal group (terminal) |
| `y_hap_isogg_freq.tsv` | Y haplogroup (ISOGG) frequency table |
| `y_hap_isogg_freq_haplists.tsv` | Individual haplogroup members per basal group (ISOGG) |
| `mt_hap_freq.tsv` | mtDNA haplogroup frequency table |
| `mt_hap_freq_haplists.tsv` | Individual haplogroup members per basal group (mtDNA) |

All frequency tables include: Ancient population name, sex, mean age (BP), country, mean latitude/longitude, per-basal-haplogroup counts, and total sample count.

### From `Haplogrpviewer.py`

| Output | Description |
|--------|-------------|
| Interactive map | Geographic locations of ancient populations with clustered markers |
| Pie chart | Basal haplogroup composition for a selected population |
| Subhaplogroup table | Tabular listing of all individual haplogroups per basal group |
| Sunburst diagram | Hierarchical visualization of haplogroup substructure |

---

## Screenshots
 
### Sidebar — Dataset Selection and Filters
Select between Y haplogroup (terminal), Y haplogroup (ISOGG), or mtDNA datasets. Filter populations by age range, country, and sex.
 
![Sidebar dataset selection and filters](screenshot_sidebar.png)
 
### Geographic Map and Basal Haplogroup Pie Chart
Ancient population locations are shown as clustered markers on an interactive map. Clicking a marker or selecting a row from the table displays the basal haplogroup composition as a pie chart.
 
![Geographic map and pie chart](screenshot_map_piechart.png)
 
### Subhaplogroup Table
For the selected population, a table lists all individual haplogroups belonging to each basal haplogroup group.
 
![Subhaplogroup table](screenshot_subhap_table.png)
 
### Haplogroup Substructure — Sunburst Diagram
A sunburst diagram visualizes the hierarchical relationship between basal haplogroups and their constituent subhaplogroups.
 
![Sunburst diagram](screenshot_sunburst.png)

## People Involved

* **Name:** Karnesh Sampath (student), Chandrashekar CR (TA), Prof. Eran Elhaik
* **Date:** 2026-03-13 / 2026-03-14

---

## Notes

* Run `Extracthaplogrps.py` first to generate all required input files before launching `Haplogrpviewer.py`.
* Use the `--force` flag in `Extracthaplogrps.py` to overwrite existing output files when rerunning the analysis.
* Rows where Sex is coded as `"c"` (child) or `"U"` (unknown) are excluded from all frequency tables and visualizations.
* The viewer caches data and map objects to improve performance — avoid manually clearing Streamlit cache unless datasets have been regenerated.
* The pipeline handles missing values encoded as `".."` and coordinate inconsistencies (comma-decimal formats) commonly found in large AADR annotation files.
* Only ancient samples (Age ≠ 0) are retained; modern samples are excluded during coordinate cleaning.
* The web application takes more computational time for large dataset.