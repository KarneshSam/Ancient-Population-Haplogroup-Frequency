import streamlit as st
import pandas as pd
import folium
from folium.plugins import MarkerCluster
import plotly.express as px
from streamlit_folium import st_folium

# file paths
y_term_path     = "/home/inf-41-2025/BINP29/Popgenetics/y_hap_ter_freq.tsv"
y_iso_path      = "/home/inf-41-2025/BINP29/Popgenetics/y_hap_isogg_freq.tsv"
mt_path         = "/home/inf-41-2025/BINP29/Popgenetics/mt_hap_freq.tsv"
y_term_sub_path = "/home/inf-41-2025/BINP29/Popgenetics/y_hap_ter_freq_haplists.tsv"
y_iso_sub_path  = "/home/inf-41-2025/BINP29/Popgenetics/y_hap_isogg_freq_haplists.tsv"
mt_sub_path     = "/home/inf-41-2025/BINP29/Popgenetics/mt_hap_freq_haplists.tsv"

# load the datasets
@st.cache_data
def load_data():
    return {
        "Y haplogroup (terminal mutation)": {
            "freq": pd.read_csv(y_term_path, sep="\t"),
            "sub":  pd.read_csv(y_term_sub_path, sep="\t"),
        },
        "Y haplogroup (ISOGG format)": {
            "freq": pd.read_csv(y_iso_path, sep="\t"),
            "sub":  pd.read_csv(y_iso_sub_path, sep="\t"),
        },
        "mtDNA haplogroup": {
            "freq": pd.read_csv(mt_path, sep="\t"),
            "sub":  pd.read_csv(mt_sub_path, sep="\t"),
        },
    }

datasets = load_data()

# Dataset selection sidebar
st.sidebar.header("Dataset Selection")
dataset_name = st.sidebar.selectbox("Select Dataset", list(datasets.keys()))
df = datasets[dataset_name]["freq"]
df_sub = datasets[dataset_name]["sub"]

