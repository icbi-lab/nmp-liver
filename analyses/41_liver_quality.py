# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:conda-2022-schneeberger-liver-scanpy]
#     language: python
#     name: conda-env-conda-2022-schneeberger-liver-scanpy-py
# ---

# %% [markdown]
# # Comparison by liver quality

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
import scanpy_helpers as sh
import pandas as pd
from nxfvars import nxfvars
import matplotlib.pyplot as plt
import altair as alt
from tqdm import tqdm
import decoupler as dc
import itertools
from pathlib import Path
import numpy as np
import seaborn as sns

# %%
sc.settings.set_figure_params(figsize=(4, 4), frameon=False)

# %%
adata_path = nxfvars.get(
    "adata_path",
    "../data/results/04_cell_type_annotation/artifacts/adata_cell_types.h5ad",
)
gene_set_il_path = nxfvars.get(
    "gene_set_il_path", "../tables/gene_sets_interleukins_chemokines.xlsx"
)

# %%
adata = sc.read_h5ad(adata_path)

# %%
gene_set_il = pd.read_excel(gene_set_il_path)

# %%
adata_t0t1 = adata[adata.obs["timepoint"].isin(["T0", "T1"]), :]

# %%
pb = sh.pseudobulk.pseudobulk(
    adata_t0t1, groupby=["patient_id", "timepoint", "cell_type"], log_norm=True
)

# %%
