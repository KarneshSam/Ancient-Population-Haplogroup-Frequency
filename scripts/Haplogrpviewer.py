#!/usr/bin/env python3

"""
Haplogrpviewer.py (script's name)

Description:
      This program is an interactive web application that visualizes ancient
      population haplogroup distributions using geographic mapping and
      interactive charts.
      
      The application loads haplogroup frequency tables and subhaplogroup
      lists generated from the AADR annotation dataset. These datasets
      summarize the distribution of Y-chromosome and mitochondrial haplogroups
      across ancient populations.

      The main objective of this application is to allow users to explore
      haplogroup frequencies geographically and inspect the internal
      structure of haplogroups within ancient populations.

      Users can filter the dataset by age range, country, and sex, and then
      interactively examine populations on a map. Selecting a population
      displays haplogroup composition as pie charts and visualizes
      subhaplogroup relationships using a sunburst diagram.

      This tool integrates geographic visualization with hierarchical
      haplogroup structure, making it easier to interpret genetic patterns
      in ancient populations.

      It also ensures robustness through error handling for missing files,
      invalid inputs, and empty datasets, providing informative messages to guide users
      and ensuring that user interactions do not cause runtime failures.
      
List of Functions:
      1. parse_arguments – Parses command-line arguments specifying
                           input dataset paths for haplogroup frequency
                           tables and subhaplogroup lists.

      2. check_file – Validates that required input files exist before
                      loading them.

      3. load_data – Loads haplogroup frequency tables and subhaplogroup
                     list files using pandas and caches them to improve
                     application performance.

      4. filter_data – Applies user-selected filters such as age range,
                       country, and sex to the frequency dataset.

      5. build_base_map – Constructs the base geographic map using
                          Folium and adds clustered markers for each
                          ancient population.

      6. Session State Logic – Maintains application state
                               including selected population,
                               map centering, and dataset changes.

List of non-standard modules:
      1. streamlit – Used to build the interactive web application
                     and manage the user interface.

      2. pandas – Used for reading datasets and performing dataframe 
                  manipulation.

      3. folium – Used for geographic visualization of ancient 
                  population locations.

      4. folium.plugins.MarkerCluster – Used to cluster nearby population 
                                        markers for improved map performance.   

      5. plotly.express – Used to create pie charts and sunburst diagrams for
                          haplogroup visualization.

      6. streamlit_folium – Used to integrate Folium maps directly into the 
                            Streamlit interface.

Procedure:
      1. The program parses command-line arguments that specify input dataset 
         paths for haplogroup frequency tables and subhaplogroup lists.
      2. Input files are validated to ensure they exist before
         loading begins.
      3. Haplogroup frequency tables and corresponding subhaplogroup list files 
         are loaded into pandas dataframes and cached to improve performance 
         during interactive exploration.
      4. Users select a dataset type from the sidebar:
            1. Y haplogroup (terminal mutation format)
            2. Y haplogroup (ISOGG format)
            3. mtDNA haplogroup
      5. Sidebar filters allow users to restrict the dataset based on:
            1. Age range
            2. Country
            3. Sex
      6. The filtered dataset is displayed in a table where users can
         select a population.
      7. The filtered populations are plotted geographically on an
         interactive map with clustered markers.
      8. Clicking a population marker on the map or selecting a row
         from the table updates the application state.
      9. The application displays:
            1. A pie chart showing basal haplogroup composition.
            2. A table listing subhaplogroups for each basal haplogroup.
            3. A sunburst diagram illustrating haplogroup substructure.
      10. The visualizations update dynamically based on user interactions 
          without reloading the entire dataset.

Input:
      1. Y haplogroup terminal frequency table (.tsv)
      2. Y haplogroup ISOGG frequency table (.tsv)
      3. mtDNA haplogroup frequency table (.tsv)
      4. Haplogroup list files corresponding to each frequency table.

Output:
      1. Interactive geographic map showing ancient population locations
      2. Basal haplogroup composition pie charts
      3. Tables listing subhaplogroups for selected populations
      4. Sunburst diagrams illustrating haplogroup hierarchy

Usage: streamlit run HaplogroupViewer.py -- \
            --y_term y_term\
            --y_iso y_iso\
            --mt mt\
            --y_term_sub y_term_sub\
            --y_iso_sub y_iso_sub\
            --mt_sub mt_sub

Get Help: streamlit run HaplogroupViewer.py -- --help


Version: 1.00
Date: 2026-03-14
Name: Karnesh Sampath

"""

# Streamlit for interactive web application
# Folium for map visualization
# Plotly for pie charts and sunburst diagrams
# argparse for command-line argument parsing
# pandas for data manipulation
# os and sys for file path handling and system interactions
import streamlit as st
import pandas as pd
import folium
from folium.plugins import MarkerCluster
import plotly.express as px
from streamlit_folium import st_folium
import argparse
import sys
import os

#######################################
# 1. ARGUMENT PARSING
#######################################

# Define command-line arguments for dataset paths
def parse_arguments():
    parser = argparse.ArgumentParser(description="Ancient Haplogroup Viewer")
    parser.add_argument("--y_term", required=True, 
                        help="Y haplogroup terminal mutation frequency")
    parser.add_argument("--y_iso", required=True,
                        help="Y haplogroup ISOGG frequency")
    parser.add_argument("--mt", required=True,
                        help="mtDNA haplogroup frequency")
    parser.add_argument("--y_term_sub", required=True,
                        help="Y terminal mutation subhaplogroup list")
    parser.add_argument("--y_iso_sub", required=True,
                        help="Y ISOGG subhaplogroup list")
    parser.add_argument("--mt_sub", required=True,
                        help="mtDNA subhaplogroup list")
    # Parse known args to allow Streamlit's own arguments to pass through
    args, _ = parser.parse_known_args()
    return args
# Parse arguments
args = parse_arguments()

# Store dataset paths in variables
# for Frequency data
Y_TERM_PATH = args.y_term
Y_ISO_PATH = args.y_iso
MT_PATH = args.mt

# for Subhaplogroup list data
Y_TERM_SUB = args.y_term_sub
Y_ISO_SUB = args.y_iso_sub
MT_SUB = args.mt_sub

#######################################
# 2. CHECK FILES EXISTENCE
#######################################

# Check if all provided file paths exist, if not, print an error and exit
def check_files(paths):
    if not os.path.exists(paths):
        st.error(f"⚠️ Error: File not found - {paths}")  
        st.stop()

check_files(Y_TERM_PATH)
check_files(Y_ISO_PATH)
check_files(MT_PATH)
check_files(Y_TERM_SUB)
check_files(Y_ISO_SUB)
check_files(MT_SUB)

#######################################
# 3. PAGE SETUP
#######################################

# Streamlit page configuration
st.set_page_config(layout="wide")
st.title("🧬 Ancient Population Haplogroup Explorer")
st.markdown("Explore ancient haplogroup distributions interactively")
#######################################
# 4. LOAD DATA
#######################################

# Load all datasets at once and cache them to avoid reloading on every interaction
# If there is an error loading any dataset, display an error message and stop the app
@st.cache_data
def load_data():
    try:
      datasets = {
        "Y haplogroup (terminal mutation)": {
            "freq": pd.read_csv(Y_TERM_PATH, sep="\t"),
            "sub":  pd.read_csv(Y_TERM_SUB, sep="\t")},
        "Y haplogroup (ISOGG format)": {
            "freq": pd.read_csv(Y_ISO_PATH, sep="\t"),
            "sub":  pd.read_csv(Y_ISO_SUB, sep="\t")},
        "mtDNA haplogroup": {
            "freq": pd.read_csv(MT_PATH, sep="\t"),
            "sub":  pd.read_csv(MT_SUB, sep="\t")}}
      return datasets
    
    except Exception as e:
        st.error(f"⚠️ Error loading datasets: {e}")
        st.stop()

# Load datasets
datasets = load_data()

#######################################
# 5. FILTERING FUNCTION
#######################################

# Apply filters to the main frequency dataframe based on user input
# create a copy of the filtered dataframe to avoid SettingWithCopyWarning
@st.cache_data
def filter_data(df, age_min, age_max, countries, sexes):
    return df[
        (df.Age >= age_min) & (df.Age <= age_max) &
        (df.Country.isin(countries)) & (df.Sex.isin(sexes))
        ].copy()

#######################################
# 6. BUILD BASE MAP FUNCTION
#######################################

# Cached map building — plain Marker objects, cached so only rebuilds on filter change
# markers added outside cache to reflect live selection state without full map rebuild
@st.cache_data
def build_base_map(df_subset, center_lat, center_long, zoom):
    m = folium.Map(location=[center_lat, center_long], zoom_start=zoom, tiles=None)
    folium.TileLayer(
        tiles="https://cartodb-basemaps-a.global.ssl.fastly.net/light_all/{z}/{x}/{y}{r}.png",
        attr="©OpenStreetMap, ©CartoDB",
        control=False,
    ).add_to(m)
    # Use MarkerCluster to group nearby markers and improve performance
    mc = MarkerCluster(options={"disableClusteringAtZoom": 10}).add_to(m)
    
    # Add markers for each row in the subset
    for _, r in df_subset.iterrows():
        marker_key = f"{r['Ancient pop']}|{r['Sex']}"
        folium.Marker(
            location=[r["Lat"], r["Long"]],
            popup=folium.Popup(marker_key, parse_html=False),
            tooltip=f"{r['Ancient pop']} | {r['Sex']} | n={r['total']}",
            icon=folium.Icon(color="blue", icon="info-sign"),
        ).add_to(mc)

    return m

#######################################
# 7. SESSION STATE AND INTERACTION LOGIC
#######################################

# Initialize session state variables to track map centering and selected population
# fly_to: (lat, long) to control map centering
# clicked_pop, clicked_sex: to track which population is selected
for key, default in [
    ("fly_to", None),
    ("clicked_pop", None),
    ("clicked_sex", None),
    ("prev_dataset", None),
    ("table_key", 0)]:
    if key not in st.session_state:
        st.session_state[key] = default

#######################################
# 8. SIDEBAR SELECTION OF DATASET 
#######################################

# Dataset selection sidebar
st.sidebar.header("📂 Dataset Selection")
st.sidebar.markdown("Choose which haplogroup dataset to explore.")
dataset_name = st.sidebar.selectbox("Select Dataset", list(datasets.keys()))

# Reset state when dataset changes
if st.session_state.prev_dataset != dataset_name:
    st.session_state.clicked_pop = None
    st.session_state.clicked_sex = None
    st.session_state.fly_to = None
    st.session_state.table_key += 1
    st.session_state.prev_dataset = dataset_name

# Get the selected dataset's frequency and subhaplogroup data
df = datasets[dataset_name]["freq"]
df_sub = datasets[dataset_name]["sub"]

#######################################
# 9. SIDEBAR FILTERS 
#######################################

# Filters sidebar
st.sidebar.header("⚙️ Filters")
st.sidebar.markdown("Filter populations by age, country, or sex.")
# create a age slider with min and max from the dataset, default to full range
age_range = st.sidebar.slider("⏳ Age Range", int(df.Age.min()), int(df.Age.max()), (int(df.Age.min()), int(df.Age.max())))
# create a multiselect for country and sex, default to all options
country_filter = st.sidebar.multiselect("📌 Country", df.Country.unique(), default=df.Country.unique())
sex_filter = st.sidebar.multiselect("⚥ Sex", df.Sex.unique(), default=df.Sex.unique())

# Apply filters to the main dataframe and 
# cache the result to avoid re-filtering on every interaction
filtered = filter_data(
    df,
    age_range[0], age_range[1],
    tuple(sorted(country_filter)),
    tuple(sorted(sex_filter)))

# If the filtered dataframe is empty, display a warning and 
# stop the app to avoid errors in the map and pie chart rendering
if filtered.empty:
    st.warning("No populations match the selected filters. Please adjust the filters to see results.")
    st.stop()

# Display filtered data
st.subheader(f"🌍 Population Table ({dataset_name})")
st.markdown("Select a population to view its geographic location and haplogroup composition.")
# user can select the population from the table
selection = st.dataframe(
    filtered.reset_index(drop=True),
    use_container_width=True,
    on_select="rerun",     # rerun the script automatically
    selection_mode="single-row")

# stores the index of the selected row 
selected_rows = selection.selection.rows
# if a row is selected, update the fly_to state to center the map on that population 
# and update the clicked population and clicked sex state variables to reflect the selection
if selected_rows:
    sel_row = filtered.iloc[selected_rows[0]]
    new_fly_to = (float(sel_row["Lat"]), float(sel_row["Long"]))
    if st.session_state.fly_to != new_fly_to:
        st.session_state.fly_to = new_fly_to
        st.session_state.clicked_pop = sel_row["Ancient pop"]
        st.session_state.clicked_sex = sel_row["Sex"]

# side by side columns for map and pie chart
col1, col2 = st.columns([2,1])

#######################################
# 10. MAP AND INTERACTION LOGIC
#######################################

# Map visualization in the first column — cached base map with dynamic marker for selection
with col1:
    st.subheader("🗺️ Geographic Distribution of Ancient Populations")
    st.markdown("View the geographic locations of ancient populations ")
    # Determine map center and zoom level based on selection or filtered data
    if st.session_state.fly_to:
        center_lat, center_long = st.session_state.fly_to
        zoom = 7
    elif not filtered.empty:
     center_lat = filtered["Lat"].mean()
     center_long = filtered["Long"].mean()
     zoom = 5
    else:
        center_lat, center_long, zoom = 20.0, 0.0, 2
    # Build base map — cached, won't rebuild unless filtered data changes
    m = build_base_map(filtered, center_lat, center_long, zoom)
    
    # Red marker added outside cache — reflects live selection state
    if st.session_state.fly_to and selected_rows:
        sel_row    = filtered.iloc[selected_rows[0]]
        marker_key = f"{sel_row['Ancient pop']}|{sel_row['Sex']}"
        folium.Marker(
            location=list(st.session_state.fly_to),
            popup=folium.Popup(marker_key, parse_html=False),
            tooltip=f"{sel_row['Ancient pop']} | {sel_row['Sex']} | n={sel_row['total']}",
            icon=folium.Icon(color="red", icon="map-marker"),
        ).add_to(m)
    # Folium map in Streamlit with a unique key to preserve state across interactions    
    map_key = (
        f"map_{dataset_name}_{age_range[0]}_{age_range[1]}_"
        f"{','.join(sorted(country_filter))}_{','.join(sorted(sex_filter))}")
    map_data = st_folium(m, width=800, height=600,
                         returned_objects=["last_object_clicked_popup"],  # capture popup clicks
                         key=map_key)  

    # Map marker click — parse popup, update state, clear table selection
    raw_popup = map_data.get("last_object_clicked_popup") if map_data else None
    if raw_popup and "|" in raw_popup:
        # Strip any HTML wrapping that folium might add
        clean = raw_popup.strip().replace("<br>", "").replace("\n", "").strip()
        if "|" in clean:
            p, s = clean.split("|", 1)
            p = p.strip()   # population
            s = s.strip()   # sex
            # Only update state if the clicked marker is different from the current selection to avoid unnecessary reruns
            if (p != st.session_state.clicked_pop or s != st.session_state.clicked_sex
                    or st.session_state.fly_to is not None):
                st.session_state.clicked_pop = p
                st.session_state.clicked_sex = s
                st.session_state.fly_to = None
                st.session_state.table_key  += 1
                st.rerun()

#######################################
# 11. PIE CHART POPULATION COMPOSITION 
#######################################

# Pie chart of basal haplogroup composition for the selected population in the second column
with col2:
    st.subheader("📊 Basal Haplogroup Composition")
    st.markdown("🔵 Click a marker on the map to see the haplogroup composition.")
    clicked_pop = st.session_state.clicked_pop
    clicked_sex = st.session_state.clicked_sex
    hap = None
    
    if clicked_pop and clicked_sex:
        # select row from filtered dataframe based on clicked population and sex
        row = df[(df["Ancient pop"] == clicked_pop) & (df["Sex"] == clicked_sex)]
        # if the row exists, get the haplogroup frequencies and create a pie chart
        if not row.empty:
             row = row.iloc[0]
             # get haplogroup frequencies for this population
             hap_cols = [c for c in df.columns if c not in ["Ancient pop","Country","Age","Lat","Long","Sex","total"]]
             hap = row[hap_cols]
             hap = hap[hap > 0]  # remove zero haplogroups
             # create pie chart using plotly express
             # only create the chart if there are haplogroups to display, otherwise show a warning
             if not hap.empty:
                fig = px.pie(values=hap, title=f"{clicked_pop} | {clicked_sex}", names=hap.index)
                
                # customize the pie chart appearance
                # show both percentage and label on the chart, and 
                # customize the hover template to show label, value, and percentage
                fig.update_traces(textinfo='percent+label',
                        hovertemplate='%{label}: %{value}<br>%{percent}',
                        textfont=dict(size=12, color="black"),
                        domain=dict(x=[0,1], y=[0,1]))
                fig.update_layout(width=350, height=350,
                        legend=dict(title=dict(text="Basal Haplogroups", font=dict(size=12, color="black")),
                        font=dict(size=10, color="black"),orientation="v",
                        yanchor="middle", y=0.5, xanchor="left", x=1.05),
                        title_font=dict(size=12, color="black", family="Arial Black"),
                        paper_bgcolor="white", plot_bgcolor="white")
                st.plotly_chart(fig)
             else:
                st.warning("⚠️ No haplogroup data for this population.")
        else:
            st.warning("⚠️ Population not found in dataset.")
    else:
        st.info("Select a population from the table or map to see its haplogroup composition.")

#######################################
# 12. SUBHAPLOGROUP TABLE AND SUNBURST DIAGRAM
#######################################

# If a population is selected, display a table of subhaplogroups and 
# a sunburst diagram showing the relationship between basal haplogroups and their subhaplogroups
if clicked_pop and hap is not None and not hap.empty:
    st.markdown("---")
    st.markdown(f"### 🧑‍🤝‍🧑 Population: {clicked_pop} ({clicked_sex})")
    st.markdown("#### 📋 Subhaplogroup lists")
    st.markdown("Tabular view of all subhaplogroups for the selected population.")
    
    # select row from subhap df
    sub_row = df_sub[(df_sub["Ancient pop"] == clicked_pop) & (df_sub["Sex"] == clicked_sex)]
    # if the row exists, create a table of basal haplogroups and their subhaplogroups, 
    # and a sunburst diagram showing the relationship between them
    if not sub_row.empty:
        sub_row = sub_row.iloc[0]   # get the first (and should be only) matching row
        # create a table
        table_data = []
        sunburst_data = {"Basal": [], "Sub": [], "Count": []}
        
        # iterate over basal haplogroups and get their subhaplogroups from the subhap dataframe
        for basal in hap.index:
            sub_text = sub_row.get(basal, "")
            # if there is no subhaplogroup data for this basal haplogroup, skip it
            if pd.isna(sub_text) or str(sub_text).strip() == "":
                continue
            # split the subhaplogroup text by comma and strip whitespace, 
            # then add to the table data and sunburst data
            sub_list = [s.strip() for s in str(sub_text).split(",")]
            table_data.append({"Basal Haplogroup": basal, "Subhaplogroups": ", ".join(sub_list)})

            # sunburst values
            count = hap[basal]
            val = count / len(sub_list)
            for s in sub_list:
                sunburst_data["Basal"].append(basal)
                sunburst_data["Sub"].append(s)
                sunburst_data["Count"].append(val)

        # display table - Subhaplogroup lists (right-aligned)
        if table_data:
            st.table(pd.DataFrame(table_data).reset_index(drop=True))
    
        # sunburst plot (left-aligned)
        if sunburst_data["Basal"]:
            st.markdown("#### 🌳 Haplogroup Substructure")
            st.markdown("Explore subhaplogroups and their hierarchy using the sunburst diagram.")
            col_sb1, col_sb2 = st.columns([1,2])
            
            with col_sb1:
                sunburst_df = pd.DataFrame(sunburst_data)
                # create sunburst diagram using plotly express, 
                # with basal haplogroups as the first level and 
                # subhaplogroups as the second level, 
                # and the count as the value
                sunburst_fig = px.sunburst(
                    sunburst_df,
                    path=['Basal','Sub'],
                    values='Count', title=f"{clicked_pop} | {clicked_sex}",
                    color='Basal')
                # customize the sunburst diagram appearance
                sunburst_fig.update_traces(textinfo='label+percent parent')
                sunburst_fig.update_layout(
                    width=450,
                    height=450, title_font=dict(size=12, color="black", family="Arial Black"),
                    paper_bgcolor="white",
                    plot_bgcolor="white")
                st.plotly_chart(sunburst_fig)
    else:
        st.info("No subhaplogroup data available for this population.")
