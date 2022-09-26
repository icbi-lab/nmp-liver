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

# %%
import scanpy as sc

# %%
adata = sc.read_h5ad("../data/results/02_scvi/artifacts/adata_scvi_doublet_filtered.h5ad")

# %%
adata.obs.groupby(["patient_id", "timepoint"]).size().reset_index(name="n").set_index(["patient_id", "timepoint"])

# %%
from tqdm.notebook import tqdm

# %%
tqdm(range(10000000))

# %%
