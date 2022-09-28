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
# # T0 vs T1
#
# functional analysis

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
from pathlib import Path
import seaborn as sns

# %%
sc.settings.set_figure_params(figsize=(4, 4), frameon=False)

# %%
adata_path = nxfvars.get(
    "adata_path",
    "../data/results/04_cell_type_annotation/artifacts/adata_cell_types.h5ad",
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/nmp-temp/")
de_res_dir = nxfvars.get("de_res_dir", "../data/results/21_deseq/DESEQ_T0_T1/")
dorothea_net = nxfvars.get("dorothea", "../tables/dorothea_human_AB_2022-09-28.csv")

# %% [markdown]
# ## Load data

# %%
adata = sc.read_h5ad(adata_path)

# %%
tfnet = pd.read_csv(dorothea_net)

# %%
adata_t0t1 = adata[adata.obs["timepoint"].isin(["T0", "T1"]), :]

# %%
de_res = {}
for f in Path(de_res_dir).glob("**/*.tsv"):
    ct = f.name.replace("_DESeq2_result.tsv", "").split("_")[-1]
    de_res[ct] = pd.read_csv(f, sep="\t")

# %%
pb = dc.get_pseudobulk(
    adata_t0t1,
    sample_col="patient_id",
    groups_col="timepoint",
    min_prop=0.05,
    min_cells=10,
    min_counts=1000,
    min_smpls=3,
)

# %%
sc.pl.umap(adata, color="cell_type")

# %% [markdown]
# ## Overview

# %%
sc.pl.umap(adata_t0t1, color="cell_type")

# %%
sc.pl.umap(adata_t0t1, color="timepoint")

# %% [markdown]
# ### cell-type fractions

# %%
cell_type_fractions = (
    adata_t0t1.obs.groupby(["patient_id", "timepoint"])
    .apply(
        lambda x: x["cell_type"].value_counts(normalize=True).rename_axis("cell_type")
    )
    .reset_index(name="frac")
)

# %%
tmp_df = (
    cell_type_fractions.groupby(["timepoint", "cell_type"]).agg("mean").reset_index()
)
alt.Chart(tmp_df).encode(x="timepoint", y="frac", color="cell_type").mark_bar()

# %%
ad_ct = sc.AnnData(
    cell_type_fractions.pivot(
        index=["patient_id", "timepoint"], columns="cell_type", values="frac"
    )
)
ad_ct.obs = ad_ct.obs.reset_index()

# %%
sh.pairwise.plot_paired(
    ad_ct,
    "timepoint",
    paired_by="patient_id",
    n_cols=8,
    panel_size=(2.4, 3.5),
    ylabel="frac",
)

# %% [markdown]
# ## Transcription factors (Dorothea)

# %%
# TFs model was obtained with this function from omnipathdb.
# For reproducibility, the result was stored in the `tables` directory and is used for all analyses.
# tfnet = dc.get_dorothea(organism="human", levels=["A", "B"])

# %%
logfc_mat = pd.concat(
    [
        r.loc[:, ["gene_id", "log2FoldChange"]]
        .rename(columns={"log2FoldChange": ct})
        .set_index("gene_id")
        for ct, r in de_res.items()
    ],
    axis=1,
    sort=True
).fillna(0).T

# %%
tf_acts, tf_pvals = dc.dense_run(dc.run_mlm, mat=logfc_mat, net=tfnet, verbose=True)

# %%
sns.clustermap(tf_acts.fillna(0), center=0, cmap="coolwarm")

# %%
