from utils.clustering_utils import seqs2hamming
import pandas as pd
import streamlit as st
import logomaker


@st.cache_data
def get_clustering_res(chain='TRB'):
    cm = pd.read_csv(f'data/significant_clone_matrix_fisher_fmba_{chain}_top_500k_wo_leaks.csv').drop(
        columns=['Unnamed: 0'])
    covid_clones_beta = cm.columns[1:]
    return seqs2hamming(covid_clones_beta, viz_method='graphopt' if chain == 'TRB' else 'drl')

@st.cache_data
def get_associated_epitopes(chain='TRB'):
    df = pd.read_csv(f'figures/associations/fmba_{chain}_sign_assoc.csv').drop(columns=['Unnamed: 0'])
    df = df[df['antigen.species'] == 'SARS-CoV-2'][['antigen.epitope', 'antigen.species', 'cluster']]
    df['chain'] = chain
    return df


@st.cache_data
def get_clonotypes_data(alpha_epitopes, beta_epitopes):
    vdjdb = pd.read_csv('data/vdjdb.txt', sep='\t')
    df = pd.read_csv(
        "figures/covid_clonotypes.csv"
    ).drop(columns=['Unnamed: 0']).merge(pd.concat([alpha_epitopes, beta_epitopes]), how='left').merge(
        vdjdb[['antigen.epitope', 'mhc.a', 'mhc.b']], how='left')
    return df


def get_coloring(clustering_res, useful_cdrs):
    useful_clusters = set(clustering_res[clustering_res.cdr3.isin(useful_cdrs)].cluster)
    clustering_res_to_plot = clustering_res.copy()
    clustering_res_to_plot['cluster'] = clustering_res_to_plot['cluster'].apply(lambda x: x if x in useful_clusters else -1)
    return clustering_res_to_plot


def plot_logo(clonotypes, ax):
    mat_df = logomaker.alignment_to_matrix(clonotypes)
    logomaker.Logo(mat_df, color_scheme='skylign_protein', ax=ax)

@st.cache_data
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')
