import os
import sys

sys.path.append(os.getcwd())
import streamlit as st

from app.prepare_data import *
from app.widgets import *

def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    st.markdown(css, unsafe_allow_html=True)
    df = df.copy()

    modification_container = main.container()
    columns_to_return = ['cdr3', 'cluster', 'SARS-CoV-2 match in vdjdb', 'chain']
    alpha_tab, beta_tab = modification_container.tabs(['TCR alpha chain clustering', 'TCR beta chain clustering'])

    with modification_container:
        df, columns_to_return = dataframe_filtration(menu, df, columns_to_return)
        filtered = df[columns_to_return].drop_duplicates().reset_index(drop=True)
    with alpha_tab:
        plot_tabs(alpha_tab, filtered, alpha_clustering, alpha_epitopes, type='TRA')
    with beta_tab:
        plot_tabs(beta_tab, filtered, beta_clustering, beta_epitopes, type='TRB')
    return df[columns_to_return].drop_duplicates().reset_index(drop=True)


st.set_page_config(layout="wide", page_title="COVID-19 TCR biomarkers")
menu, mid_margin, main, margin_right = st.columns((10, 1, 20, 2))

create_side_menu()

beta_clustering = get_clustering_res('TRB')
alpha_clustering = get_clustering_res('TRA')
alpha_epitopes = get_associated_epitopes('TRA')
beta_epitopes = get_associated_epitopes('TRB')

df = get_clonotypes_data(alpha_epitopes, beta_epitopes).rename(columns={'has_covid_association': 'SARS-CoV-2 match in vdjdb'})

filtered = filter_dataframe(df)

