import os
import sys

sys.path.append(os.getcwd())
from app.params import *
import pandas as pd
import streamlit as st
from pandas.api.types import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_object_dtype,
)
from utils.clustering_utils import check_distance
from app.prepare_data import *
import plotly.express as px
import matplotlib.pyplot as plt

st.set_page_config(layout="wide")
st.title("Predicting SARS-CoV-2 exposure using T-cell repertoire sequencing and machine learning")
menu, mid_margin, main, margin_right = st.columns((10, 2, 20, 2))

github_b, zenodo_b, preprint_b = menu.columns((1, 1, 1))
github_b.markdown(
    f'<a href="{github_link}" style="{button_style}">GitHub repo</a>',
    unsafe_allow_html=True
)
zenodo_b.markdown(
    f'<a href="{zenodo_link}" style="{button_style}">Zenodo record</a>',
    unsafe_allow_html=True
)
preprint_b.markdown(
    f'<a href="{preprint_link}" style="{button_style}">Preprint</a>',
    unsafe_allow_html=True
)
menu.markdown(abstract_text, unsafe_allow_html=True)
menu.markdown(desc_text, unsafe_allow_html=True)


def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:

    df = df.copy()

    for col in df.columns:
        if is_object_dtype(df[col]):
            try:
                df[col] = pd.to_datetime(df[col])
            except Exception:
                pass

        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].dt.tz_localize(None)

    modification_container = main.container()
    columns_to_return = ['cdr3', 'cluster', 'has_covid_association', 'chain']

    with modification_container:
        to_filter_columns = ['cdr3', 'chain', 'has_covid_association','mhc.a', 'mhc.b']
        for column in to_filter_columns:
            left, mid = menu.columns((1, 15))
            left.write("↳")
            if is_categorical_dtype(df[column]) or df[column].nunique() < 10 or 'mhc' in column:
                user_cat_input = mid.multiselect(
                    f"Values for {column}",
                    df[column].unique(),
                    default=[],
                )
                if len(user_cat_input) == 0:
                    user_cat_input = list(df[column].unique())
                elif column not in columns_to_return:
                        columns_to_return.append(column)
                df = df[df[column].isin(user_cat_input)]
            elif is_numeric_dtype(df[column]):
                _min = float(df[column].min())
                _max = float(df[column].max())
                step = (_max - _min) / 100
                user_num_input = mid.slider(
                    f"Values for {column}",
                    _min,
                    _max,
                    (_min, _max),
                    step=step,
                )
                df = df[df[column].between(*user_num_input)]
            elif is_datetime64_any_dtype(df[column]):
                user_date_input = mid.date_input(
                    f"Values for {column}",
                    value=(
                        df[column].min(),
                        df[column].max(),
                    ),
                )
                if len(user_date_input) == 2:
                    user_date_input = tuple(map(pd.to_datetime, user_date_input))
                    start_date, end_date = user_date_input
                    df = df.loc[df[column].between(start_date, end_date)]
            else:
                option = mid.radio('Perform clonotype search by:', ['regex', 'Levenstein distance'], key='opt',
                                   horizontal=True)

                if option == 'regex':
                    user_text_input = mid.text_input(
                        f"Regex corresponding to {column}",
                    )
                    if user_text_input:
                        df = df[df[column].str.contains(user_text_input)]
                else:
                    user_text_input = mid.text_input(
                        f"Sequence of {column}, evaluation by Levenstein distance",
                    )
                    mismatch = mid.select_slider(
                        f"Maximum mismatches in CDR3 sequence",
                        options=[x for x in range(5)],
                        value=1
                    )
                    if user_text_input:
                        df = df[df[column].apply(lambda x: check_distance(x, user_text_input, dist=mismatch))]
    plot_clustering_res(get_coloring(alpha_clustering, df.cdr3).merge(alpha_epitopes, how='left').fillna('none'), alpha)
    plot_clustering_res(get_coloring(beta_clustering, df.cdr3).merge(beta_epitopes, how='left').fillna('none'), beta)
    return df[columns_to_return].drop_duplicates().reset_index(drop=True)


def plot_clustering_res(df, col=None):
    df['cluster'] = df['cluster'].astype(str)
    fig = px.scatter(
        df[df.cluster != '-1'],
        x="x",
        y="y",
        color="cluster",
        hover_name="cdr3",
        hover_data=["antigen.epitope", 'antigen.species'],
        width=600, height=600
    )
    fig.add_trace(px.scatter(
        df[df.cluster == '-1'],
        x="x",
        y="y",
        hover_name="cdr3",
        hover_data=["antigen.epitope", 'antigen.species'],
        color_discrete_sequence=["grey"],
    ).data[0])
    artist = main if col is None else col
    artist.plotly_chart(fig, theme="streamlit", use_container_width=True,
                        color_discrete_sequence=px.colors.qualitative.Alphabet)


beta_clustering = get_clustering_res('TRB')
alpha_clustering = get_clustering_res('TRA')
alpha_epitopes = get_associated_epitopes('TRA')
beta_epitopes = get_associated_epitopes('TRB')

df = get_clonotypes_data(alpha_epitopes, beta_epitopes)

alpha, beta = main.columns((1, 1))
alpha.subheader('TCR alpha chain clustering')
beta.subheader('TCR beta chain clustering')
filtered = filter_dataframe(df)
filtered_alpha = filtered[filtered.chain == 'TRA']
filtered_beta = filtered[filtered.chain == 'TRB']
alpha.dataframe(filtered_alpha)
beta.dataframe(filtered_beta)

if len(filtered_alpha):
    alpha_cluster = alpha.selectbox('Plot logo for alpha cluster: ',
                                    (x for x in sorted(filtered_alpha.cluster.unique())),
                                    )
    fig, ax = plt.subplots(figsize=(5, 2))
    plot_logo(list(filtered_alpha[filtered_alpha.cluster == alpha_cluster].cdr3), ax)
    ax.yaxis.set_visible(False)
    alpha.pyplot(fig)

if len(filtered_beta):
    beta_cluster = beta.selectbox('Plot logo for beta cluster: ',
                                  (x for x in sorted(filtered_beta.cluster.unique())),
                                  )
    fig, ax = plt.subplots(figsize=(5, 2))
    plot_logo(list(filtered_beta[filtered_beta.cluster == beta_cluster].cdr3), ax)
    ax.yaxis.set_visible(False)
    beta.pyplot(fig)

