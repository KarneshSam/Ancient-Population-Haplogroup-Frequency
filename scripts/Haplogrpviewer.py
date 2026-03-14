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

    args = parser.parse_args()
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

# Apply filters
# create a copy of the filtered dataframe to avoid SettingWithCopyWarning
filtered = df[(df.Age >= age_range[0]) &
              (df.Age <= age_range[1]) &
              (df.Country.isin(country_filter)) &
              (df.Sex.isin(sex_filter))].copy()

# Display filtered data
st.subheader(f"Population Table ({dataset_name})")
st.dataframe(filtered, use_container_width=True)

# side by side columns for map and pie chart
col1, col2 = st.columns([2,1])
clicked_pop = None

# map
with col1:
    st.subheader("Geographic Distribution of Ancient Populations")

    center_lat = filtered["Lat"].mean()
    center_long = filtered["Long"].mean()
    
    m = folium.Map(location=[center_lat, center_long], zoom_start=5, tiles=None)
    folium.TileLayer(
        tiles="https://cartodb-basemaps-a.global.ssl.fastly.net/light_all/{z}/{x}/{y}{r}.png",
        attr="©OpenStreetMap, ©CartoDB", control=False).add_to(m)
    
    marker_cluster = MarkerCluster().add_to(m)
    for idx, r in filtered.iterrows():
        # unique key combining population and sex
        marker_key = f"{r['Ancient pop']}|{r['Sex']}"
        tooltip_text = f"{r['Ancient pop']} | Sex: {r['Sex']} | Total: {r['total']}"
        folium.Marker(location=[r["Lat"], r["Long"]],
                      popup=marker_key,
                      tooltip=tooltip_text,
                      icon=folium.Icon(color="blue", icon="info-sign")
                     ).add_to(marker_cluster)
        
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

