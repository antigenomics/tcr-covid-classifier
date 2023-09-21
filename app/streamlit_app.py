import os
import sys
working_dir = 'C:/Users/lizzka239/projects/tcr-covid-classifier'
os.chdir(working_dir)
sys.path.append(os.getcwd())
import pandas as pd
import streamlit as st
from pandas.api.types import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
    is_object_dtype,
)
def check_distance(cdr1, cdr2, dist=1):
    if len(cdr1) != len(cdr2):
        return False
    found_diff_count = 0
    for c1, c2 in zip(cdr1, cdr2):
        if c1 != c2 and found_diff_count < dist:
            found_diff_count += 1
        elif c1 != c2 and found_diff_count == dist:
            return False
    return True
st.title("Predicting SARS-CoV-2 exposure using T-cell repertoire sequencing and machine learning")
st.write(
    """This page accomodated the COVID-19 associated clonotypes found in the research [add link here!]. You can 
    search for the clonotypes which are similar to the desired one and visualize the closest clonotypes.
    """
)


def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a UI on top of a dataframe to let viewers filter columns

    Args:
        df (pd.DataFrame): Original dataframe

    Returns:
        pd.DataFrame: Filtered dataframe
    """
    modify = st.checkbox("Add filters")

    if not modify:
        return df

    df = df.copy()

    # Try to convert datetimes into a standard format (datetime, no timezone)
    for col in df.columns:
        if is_object_dtype(df[col]):
            try:
                df[col] = pd.to_datetime(df[col])
            except Exception:
                pass

        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].dt.tz_localize(None)

    modification_container = st.container()

    with modification_container:
        to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
        for column in to_filter_columns:
            left, mid, right = st.columns((1, 15, 5))
            left.write("↳")
            # Treat columns with < 10 unique values as categorical
            if is_categorical_dtype(df[column]) or df[column].nunique() < 10:
                user_cat_input = mid.multiselect(
                    f"Values for {column}",
                    df[column].unique(),
                    default=list(df[column].unique()),
                )
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
                user_text_input = mid.text_input(
                    f"Substring or regex in {column}",
                )
                mismatch = right.select_slider(
                    f"Allow X mismatches in CDR sequence",
                    options=[x for x in range(5)],
                    value=1
                )
                if user_text_input:
                    df = df[df[column].apply(lambda x: check_distance(x, user_text_input, dist=mismatch))]
                    # df = df[df[column].str.contains(user_text_input)]

    return df


df = pd.read_csv(
    "figures/covid_clonotypes.csv"
).drop(columns=['Unnamed: 0'])
st.dataframe(filter_dataframe(df))