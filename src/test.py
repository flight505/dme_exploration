# %% [markdown]
"""

"""
# %%
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.io as pio

# %%


def load_data(path: str = "DB_MTX_USA.xlsx"):
    return pd.read_excel(path).sort_values(["Patient_id", "Sample_time"])


# %%
df = load_data()
df
# %%
# Load and sort data

dose_mtx_df = df.loc[
    df["Sample_type"] == "Dose_MTX", ["Patient_id", "Sample_time", "Result"]
]
dose_mtx_df

# %%
def _to_number(x):
    """Remove the '<0.0x' strings from result..."""
    try:
        return float(x)
    except:
        return 0.0


# %%
dose_mtx_df["Result"] = dose_mtx_df["Result"].apply(_to_number)
dose_mtx_df
# %%
