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
# # Overview plots

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
import scanpy_helpers as sh
import pandas as pd
import numpy as np
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
cell_type_marker_path = nxfvars.get(
    "marker_genes_path", "../tables/cell_type_markers.csv"
)

# %% [markdown]
# ## Load data

# %%
cell_type_markers = pd.read_csv(cell_type_marker_path)

# %%
ah = sh.annotation.AnnotationHelper(markers=cell_type_markers.sort_values("cell_type"))

# %%
adata = sc.read_h5ad(adata_path)

# %%
pb_cell_type_coarse = dc.get_pseudobulk(
    adata,
    sample_col="patient_id",
    groups_col="cell_type_coarse",
    min_prop=0.05,
    min_cells=10,
    min_counts=1000,
    min_smpls=3,
)
sc.pp.normalize_total(pb_cell_type_coarse, target_sum=1e6)
sc.pp.log1p(pb_cell_type_coarse)

# %% [markdown]
# ## UMAP plots

# %%
for col in ["cell_type", "cell_type_coarse"]:
    with plt.rc_context({"figure.figsize": (6, 6), "figure.dpi": 1200}):
        sh.colors.set_scale_anndata(adata, col)
        fig = sc.pl.umap(adata, color=col, return_fig=True, size=10)
        fig.savefig(f"{artifact_dir}/umap_{col}.pdf", bbox_inches="tight")

# %% [markdown]
# ## Cell-type markers

# %%
fig = ah.plot_dotplot(adata, groupby="cell_type", return_fig=True)
fig.savefig(f"{artifact_dir}/dotplot_cell_types.pdf", bbox_inches="tight")

# %%
# %%capture
for gene in [
    "PTPRC",
    "EPCAM",
    "SPARC",
    "FCGR3B",
    "CD68",
    "CD3E",
    "CD4",
    "CD8A",
    "FOXP3",
    "NKG7",
    "FLT3",
    "CD24",
    "CD79A",
    "JCHAIN",
    "ALB",
    "KRT18",
]:
    with plt.rc_context({"figure.figsize": (3, 3), "figure.dpi": 600}):
        fig = sc.pl.umap(adata, color=gene, return_fig=True)
        fig.savefig(f"{artifact_dir}/umap_{gene}.pdf", bbox_inches="tight")

# %% [markdown]
# ## Cell-type fractions

# %%
per_patient = (
    adata.obs.groupby(["patient_id", "cell_type"]).size().reset_index(name="n")
)

# %%
ch = (
    alt.Chart(per_patient)
    .mark_bar()
    .encode(
        x=alt.X("n", stack="normalize"),
        y="patient_id",
        color=alt.Color("cell_type", scale=sh.colors.altair_scale("cell_type")),
    )
)
ch.save(f"{artifact_dir}/cell_type_barchart_per_patient.svg")
ch

# %%
per_sample = (
    adata.obs.groupby(["patient_id", "sample_id", "timepoint", "cell_type_coarse"])
    .size()
    .reset_index(name="n")
)

# %%
per_sample.to_csv(f"{artifact_dir}/cell_type_coarse_counts_per_sample.csv")

# %% [markdown]
# ## Transcript counts

# %%
with plt.rc_context({"figure.figsize": (6, 6), "figure.dpi": 1200}):
    fig = sc.pl.umap(
        adata, color="total_counts", vmax=20000, cmap="Reds", size=15, return_fig=True
    )
    fig.savefig(f"{artifact_dir}/umap_total_counts.pdf", bbox_inches="tight")

# %%
mean_counts = (
    adata.obs.groupby(["patient_id", "cell_type_coarse"])
    .agg(total_counts=pd.NamedAgg("total_counts", "mean"))
    .reset_index()
)

# %%
order = (
    mean_counts.groupby("cell_type_coarse")
    .agg("median")
    .sort_values("total_counts")
    .index.tolist()
)

# %%
ch = (
    alt.Chart(mean_counts)
    .mark_boxplot()
    .encode(
        x=alt.X("cell_type_coarse", sort=order),
        y=alt.Y("total_counts"),
        color=alt.Color(
            "cell_type_coarse", scale=sh.colors.altair_scale("cell_type_coarse")
        ),
    )
)
ch.save(f"{artifact_dir}/transcript_counts_per_cell_type.svg")

# %% [markdown]
# ## Gene expression across cell-types

# %%
tmp_df = pb_cell_type_coarse.obs
tmp_df["CXCR2"] = np.array(pb_cell_type_coarse[:, "CXCR2"].X[:, 0])

# %%
order = (
    tmp_df.groupby("cell_type_coarse")
    .agg("median")
    .sort_values("CXCR2", ascending=False)
    .index.values
)

# %%
PROPS = {
    "boxprops": {"facecolor": "none", "edgecolor": "darkgrey"},
    "medianprops": {"color": "darkgrey"},
    "whiskerprops": {"color": "darkgrey"},
    "capprops": {"color": "darkgrey"},
}

fig, ax = plt.subplots(1, 1, figsize=(7, 4))
sns.stripplot(
    x="cell_type_coarse",
    y="CXCR2",
    hue="patient_id",
    data=tmp_df,
    ax=ax,
    order=order,
    palette=sh.colors.COLORS.patient_id,
    size=7,
    linewidth=1
)
sns.boxplot(
    x="cell_type_coarse",
    y="CXCR2",
    ax=ax,
    data=tmp_df,
    order=order,
    color="white",
    **PROPS,
    showfliers=False,
)
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
_ = plt.xticks(rotation=90)
# fig.savefig(f"{artifact_dir}/vegfa_fractions.pdf", bbox_inches="tight")

# %%
