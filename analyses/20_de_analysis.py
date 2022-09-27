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
# # Prepare for DE analysis
#
# Generate pseudobulk and export for analysis with DESeq2 script

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

# %%
sc.settings.set_figure_params(figsize=(4, 4), frameon=False)

# %%
adata_path = nxfvars.get(
    "adata_path",
    "../data/results/04_cell_type_annotation/artifacts/adata_cell_types.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/nmp-temp/")

# %%
adata = sc.read_h5ad(adata_path)

# %%
sc.pl.umap(adata, color="cell_type")

# %% [markdown]
# ## Make pseudobulk

# %%
pd.set_option("display.max_columns", None)

# %%
adata.obs


# %%
def save_pseudobulk(pb, samplesheet_filename, counts_filename):
    samplesheet = pb.obs.copy()
    samplesheet.drop(columns=["sample"], errors="ignore")
    samplesheet.index.name = "sample"
    samplesheet.reset_index(inplace=True)
    bulk_df = pb.to_df().T
    bulk_df.index.name = "gene_id"
    samplesheet.to_csv(samplesheet_filename)
    bulk_df.to_csv(counts_filename)


# %% [markdown]
# ### T0 vs T1

# %%
for ct, tmp_ad in tqdm(sh.util.split_anndata(adata[adata.obs["timepoint"].isin(["T0", "T1"]), :], "cell_type")):
    pb = dc.get_pseudobulk(
        tmp_ad,
        sample_col="patient_id",
        groups_col="timepoint",
        min_prop=0.05,
        min_cells=10,
        min_counts=1000,
        min_smpls=3
    )
    if pb.obs["timepoint"].nunique() <= 1:
        print(f"Cell type {ct} does not have enough replicates per group")
    else:
        save_pseudobulk(pb, f"{artifact_dir}/t0_vs_t1_{ct}.samplesheet.csv", f"{artifact_dir}/t0_vs_t1_{ct}.counts.csv")        

# %% [markdown]
# ### Transplanted vs. discarded and HQ vs marginal for all timepoints

# %%
for timepoint in ["T0", "T1", "T2"]: 
    for groups_col in ["LT", "ECD"]: 
        for ct, tmp_ad in tqdm(sh.util.split_anndata(adata[adata.obs["timepoint"] == timepoint, :], "cell_type")):
            pb = dc.get_pseudobulk(
                tmp_ad,
                sample_col="patient_id",
                groups_col=groups_col,
                min_prop=0.05,
                min_cells=10,
                min_counts=1000,
                min_smpls=2
            )
            if pb.obs[groups_col].nunique() <= 1:
                print(f"Cell type {ct} does not have enough replicates per group")
            else:
                basename = f"{artifact_dir}/{groups_col}_{timepoint}_{ct}"
                save_pseudobulk(pb, f"{basename}.samplesheet.csv", f"{basename}.counts.csv")   

# %%
