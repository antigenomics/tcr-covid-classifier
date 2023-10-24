import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px
import streamlit as st
from pandas.api.types import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
)

from app.params import *
from app.prepare_data import get_coloring, convert_df, plot_logo
from utils.clustering_utils import check_distance


def create_side_menu():
    with st.sidebar:
        st.title("Predicting SARS-CoV-2 exposure using T-cell repertoire sequencing and machine learning")
        github_b, zenodo_b, preprint_b = st.columns((1, 1, 1))
        github_b.markdown(
            f'<a href="{github_link}" style="{button_style}">GitHub</a>',
            unsafe_allow_html=True
        )
        zenodo_b.markdown(
            f'<a href="{zenodo_link}" style="{button_style}">Zenodo</a>',
            unsafe_allow_html=True
        )
        preprint_b.markdown(
            f'<a href="{preprint_link}" style="{button_style}">Preprint</a>',
            unsafe_allow_html=True
        )
        st.header('Summary')
        st.success(abstract_text, )
        st.header('Database description')
        st.success(desc_text)


def dataframe_filtration(menu, df, columns_to_return):
    to_filter_columns = ['cdr3', 'chain', 'SARS-CoV-2 match in vdjdb', 'mhc.a', 'mhc.b']
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
    return df, columns_to_return


def plot_clustering_res(df, col):
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
    artist = col
    artist.plotly_chart(fig, theme="streamlit", use_container_width=True,
                        color_discrete_sequence=px.colors.qualitative.Alphabet)


def plot_tabs(tab, filtered, clustering, epitopes, type='TRA'):
    plot_tab, logo_tab, dataframe_tab = tab.tabs(['TCR similarity map', 'CDR3 motifs', 'Data table'])

    filtered_data_of_type = filtered[filtered.chain == type]
    with plot_tab:
        plot_clustering_res(get_coloring(clustering, filtered.cdr3).merge(epitopes,
                                                                    how='left').fillna('none'),
                            plot_tab)

    with dataframe_tab:
        dataframe_tab.dataframe(filtered_data_of_type)
        dataframe_tab.download_button(
            label=f"Download selected {type} clones",
            data=convert_df(filtered_data_of_type),
            file_name=f'{type}_clones_selected.csv',
            mime='text/csv',
        )

    if len(filtered_data_of_type):
        with logo_tab:
            alpha_cluster = logo_tab.selectbox(f'Plot logo for {type} cluster: ',
                                               (x for x in sorted(filtered_data_of_type.cluster.unique())),
                                               )
            fig, ax = plt.subplots(figsize=(5, 2))
            plot_logo(list(filtered_data_of_type[filtered_data_of_type.cluster == alpha_cluster].cdr3), ax)
            ax.yaxis.set_visible(False)
            logo_tab.pyplot(fig)
