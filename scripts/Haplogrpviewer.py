import streamlit as st
import pandas as pd
import folium
from folium.plugins import MarkerCluster
import plotly.express as px
from streamlit_folium import st_folium
import argparse

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

    args, _ = parser.parse_known_args()
    return args

args = parse_arguments()

Y_TERM_PATH = args.y_term
Y_ISO_PATH  = args.y_iso
MT_PATH     = args.mt

Y_TERM_SUB  = args.y_term_sub
Y_ISO_SUB   = args.y_iso_sub
MT_SUB      = args.mt_sub

# Page config
st.set_page_config(layout="wide")
st.title("Ancient Population Haplogroup Explorer")

# load the datasets
@st.cache_data
def load_data():
    return {
        "Y haplogroup (terminal mutation)": {
            "freq": pd.read_csv(Y_TERM_PATH, sep="\t"),
            "sub":  pd.read_csv(Y_TERM_SUB, sep="\t")},
        "Y haplogroup (ISOGG format)": {
            "freq": pd.read_csv(Y_ISO_PATH, sep="\t"),
            "sub":  pd.read_csv(Y_ISO_SUB, sep="\t")},
        "mtDNA haplogroup": {
            "freq": pd.read_csv(MT_PATH, sep="\t"),
            "sub":  pd.read_csv(MT_SUB, sep="\t")},
        }

datasets = load_data()

# Apply filters
# @st.cache_data store in memory
# create a copy of the filtered dataframe to avoid SettingWithCopyWarning
@st.cache_data
def filter_data(df, age_min, age_max, countries, sexes):
    return df[
        (df.Age >= age_min) & (df.Age <= age_max) &
        (df.Country.isin(countries)) & (df.Sex.isin(sexes))
    ].copy()

# Cached map building — plain Marker objects, cached so only rebuilds on filter change
@st.cache_data
def build_base_map(df_subset, center_lat, center_long, zoom):
    m = folium.Map(location=[center_lat, center_long], zoom_start=zoom, tiles=None)
    folium.TileLayer(
        tiles="https://cartodb-basemaps-a.global.ssl.fastly.net/light_all/{z}/{x}/{y}{r}.png",
        attr="©OpenStreetMap, ©CartoDB",
        control=False,
    ).add_to(m)

    mc = MarkerCluster(options={"disableClusteringAtZoom": 10}).add_to(m)
    
    for _, r in df_subset.iterrows():
        marker_key = f"{r['Ancient pop']}|{r['Sex']}"
        folium.Marker(
            location=[r["Lat"], r["Long"]],
            popup=folium.Popup(marker_key, parse_html=False),
            tooltip=f"{r['Ancient pop']} | {r['Sex']} | n={r['total']}",
            icon=folium.Icon(color="blue", icon="info-sign"),
        ).add_to(mc)

    return m

# need a storage
# helpful when selecting a population from table
for key, default in [
    ("fly_to", None),
    ("clicked_pop", None),
    ("clicked_sex", None),
    ("prev_dataset", None),

]:
    if key not in st.session_state:
        st.session_state[key] = default

# Dataset selection sidebar
st.sidebar.header("Dataset Selection")
dataset_name = st.sidebar.selectbox("Select Dataset", list(datasets.keys()))
df = datasets[dataset_name]["freq"]
df_sub = datasets[dataset_name]["sub"]

# Filters sidebar
st.sidebar.header("Filters")
age_range = st.sidebar.slider("Age Range", int(df.Age.min()), int(df.Age.max()), (int(df.Age.min()), int(df.Age.max())))
country_filter = st.sidebar.multiselect("Country", df.Country.unique(), default=df.Country.unique())
sex_filter = st.sidebar.multiselect("Sex", df.Sex.unique(), default=df.Sex.unique())

filtered = filter_data(
    df,
    age_range[0], age_range[1],
    tuple(sorted(country_filter)),
    tuple(sorted(sex_filter))
)

# Display filtered data
st.subheader(f"Population Table ({dataset_name})")

# user can select the population from the table
selection = st.dataframe(
    filtered.reset_index(drop=True),
    use_container_width=True,
    on_select="rerun",     # rerun the script automatically
    selection_mode="single-row"
)
# stores the index of the selected row 
selected_rows = selection.selection.rows
if selected_rows:
    sel_row = filtered.iloc[selected_rows[0]]
    new_fly_to = (float(sel_row["Lat"]), float(sel_row["Long"]))
    if st.session_state.fly_to != new_fly_to:
        st.session_state.fly_to = new_fly_to
        st.session_state.clicked_pop = sel_row["Ancient pop"]
        st.session_state.clicked_sex = sel_row["Sex"]

# side by side columns for map and pie chart
col1, col2 = st.columns([2,1])

# map
with col1:
    st.subheader("Geographic Distribution of Ancient Populations")
    
    if st.session_state.fly_to:
        center_lat, center_long = st.session_state.fly_to
        zoom = 7
    elif not filtered.empty:
     center_lat = filtered["Lat"].mean()
     center_long = filtered["Long"].mean()
     zoom = 5
    # Build base map — cached, won't rebuild unless filtered data changes
    m = build_base_map(filtered, center_lat, center_long, zoom)
    
    # Red marker added outside cache — reflects live selection state
    if st.session_state.fly_to and selected_rows:
        sel_row    = filtered.iloc[selected_rows[0]]
        marker_key = f"{sel_row['Ancient pop']}|{sel_row['Sex']}"
        folium.Marker(
            location=list(st.session_state.fly_to),
            popup=folium.Popup(marker_key, parse_html=False),
            tooltip=f"{sel_row['Ancient pop']} | {sel_row['Sex']} | n={r['total']}",
            icon=folium.Icon(color="red", icon="map-marker"),
        ).add_to(m)

    map_data = st_folium(m, width=800, height=600)

# pie chart
with col2:
    st.subheader("Basal Haplogroup Distribution")
    clicked_pop = None
    clicked_sex = None
    
    if map_data and map_data.get("last_object_clicked"):
        clicked_key = map_data["last_object_clicked"]
        clicked_pop, clicked_sex = clicked_key.split("|")
    
    if clicked_pop:
        # select row from main df
        row = df[(df["Ancient pop"] == clicked_pop) & (df["Sex"] == clicked_sex)].iloc[0]
        # get haplogroup frequencies for this population
        hap_cols = [c for c in df.columns if c not in ["Ancient pop","Country","Age","Lat","Long","Sex","total"]]
        hap = row[hap_cols]
        hap = hap[hap > 0]  # remove zero haplogroups
    
    if not hap.empty:
            fig = px.pie(values=hap, title=clicked_key, names=hap.index)
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
        st.info("Click a marker on the map to see the pie chart.")

# subhaplogroup table
if clicked_pop:
    st.markdown("---")
    st.markdown(f"### Population: {clicked_pop} ({clicked_sex})")
    st.markdown("#### Subhaplogroup Breakdown")
    
    # select row from subhap df
    sub_row = df_sub[(df_sub["Ancient pop"] == clicked_pop) & (df_sub["Sex"] == clicked_sex)].iloc[0]
    # create a table
    table_data = []
    sunburst_data = {"Basal": [], "Sub": [], "Count": []}

    for basal in hap.index:
            sub_text = sub_row.get(basal, "")
            if pd.isna(sub_text) or str(sub_text).strip() == "":
                continue
            sub_list = [s.strip() for s in str(sub_text).split(",")]
            table_data.append({"Basal Haplogroup": basal, "Subhaplogroups": ", ".join(sub_list)})

            # sunburst values
            count = hap[basal]
            val = count / len(sub_list)
            for s in sub_list:
                sunburst_data["Basal"].append(basal)
                sunburst_data["Sub"].append(s)
                sunburst_data["Count"].append(val)

    # display table
    if table_data:
            st.table(pd.DataFrame(table_data))
    
    # sunburst (left-aligned)
    if sunburst_data["Basal"]:
            st.markdown("#### Haplogroup Substructure")
            col_sb1, col_sb2 = st.columns([1,2])

            with col_sb1:
                sunburst_df = pd.DataFrame(sunburst_data)
                sunburst_fig = px.sunburst(
                    sunburst_df,
                    path=['Basal','Sub'],
                    values='Count',
                    color='Basal'
                )
                sunburst_fig.update_traces(textinfo='label+percent parent')
                sunburst_fig.update_layout(
                    width=450,
                    height=450,
                    paper_bgcolor="white",
                    plot_bgcolor="white"
                )
                st.plotly_chart(sunburst_fig)

