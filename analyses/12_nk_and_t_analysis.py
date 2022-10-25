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
# %load_ext autoreload
# %autoreload 2

# %%
# ensure reproducibility -> set numba multithreaded mode
from nxfvars import nxfvars
from threadpoolctl import threadpool_limits
import os

cpus = nxfvars.get("cpus", 2)
os.environ["NUMBA_NUM_THREADS"] = str(cpus)
threadpool_limits(cpus)

# %%
import scanpy as sc
import scanpy_helpers as sh
import pandas as pd
import matplotlib.pyplot as plt
import altair as alt
from tqdm import tqdm

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
# ## Marker dotplot

# %%
t_markers = [
    "CD3E",
    "IL32",
    "TRBC2",
    "TRAC",
    "CD2",
    "ETS1",
    "RUNX3",
    "CCL5",
    "IKZF3",
    "FYN",
    "CD3D",
    "CD3G",
    "CD69",
    "FYN",
    "CD8A",
]
fig = sc.pl.dotplot(adata, groupby="cell_type_coarse", var_names=t_markers, return_fig=True, cmap="coolwarm")
fig.savefig(f"{artifact_dir}/nkt_dotplot.pdf", bbox_inches="tight")

# %% [markdown]
# ## Cluster overview

# %%
adata_nkt = adata[adata.obs["cell_type_coarse"].isin(["NK cells", "T cells"]), :].copy()

# %%
adata_nkt.shape

# %%
# Only T0 and T1
adata_nkt = adata_nkt[adata_nkt.obs["timepoint"] != "T2", :].copy()

# %%
adata_nkt.shape

# %%
sc.pp.neighbors(adata_nkt, use_rep="X_scVI", n_neighbors=40)

# %%
sc.tl.umap(adata_nkt)

# %%
for col in ["cell_type", "cell_type_coarse", "patient_id", "timepoint", "sample_id"]:
    if col != "sample_id":
        sh.colors.set_scale_anndata(adata_nkt, col)
    legend_params = (
        dict(legend_loc="on data", legend_fontoutline=2)
        if col == "cell_type" or col == "cell_type_coarse"
        else {}
    )
    with plt.rc_context({"figure.figsize": (6, 6)}):
        fig = sc.pl.umap(
            adata_nkt, color=col, return_fig=True, size=20, **legend_params
        )
        fig.savefig(
            f"{artifact_dir}/umap_nkt_cluster_overview_{col}.pdf",
            bbox_inches="tight",
            dpi=1200,
        )

# %%
for gene in ["KLRD1", "CD3E"]:
    with plt.rc_context({"figure.figsize": (6, 6)}):
        fig = sc.pl.umap(
            adata_nkt, color=gene, return_fig=True, size=20, cmap="inferno"
        )
        fig.savefig(
            f"{artifact_dir}/umap_nkt_cluster_overview_{col}_{gene}.pdf",
            bbox_inches="tight",
            dpi=1200,
        )

# %%
adata_nkt.obs["cell_type"].value_counts()

# %%
tmp_df = (
    adata_nkt.obs.groupby(["patient_id", "cell_type"])
    .size()
    .reset_index(name="n_cells")
)
heatmp = (
    alt.Chart(tmp_df)
    .mark_rect()
    .encode(
        x="cell_type",
        y=alt.Y("patient_id"),
        color=alt.Color("n_cells", scale=alt.Scale(scheme="inferno", reverse=True)),
    )
)
txt = (
    alt.Chart(tmp_df)
    .mark_text()
    .encode(
        x="cell_type",
        y=alt.Y("patient_id"),
        text="n_cells",
        color=alt.condition(
            alt.datum.n_cells < 500, alt.value("black"), alt.value("white")
        ),
    )
)
ch = (heatmp + txt).properties(width=200)
ch.save(f"{artifact_dir}/nk_and_t_cluster_count_per_patient.svg")
ch.display()

# %%
tmp_df = (
    adata_nkt.obs.groupby(["patient_id", "cell_type", "timepoint"])
    .size()
    .reset_index(name="n_cells")
    .groupby(["cell_type", "timepoint"])
    .agg("mean")
    .reset_index()
)

# %%
ch = (
    alt.Chart(tmp_df)
    .encode(
        color=alt.Color("cell_type", scale=sh.colors.altair_scale("cell_type")),
        y="timepoint",
        x=alt.X("n_cells", stack="normalize"),
    )
    .mark_bar()
)
ch.save(f"{artifact_dir}/timepoints_by_cell_type.svg")
ch.display()

# %%
ch = (
    alt.Chart(tmp_df)
    .encode(
        color=alt.Color("timepoint", scale=sh.colors.altair_scale("timepoint")),
        y=alt.Y("n_cells", stack="normalize"),
        x="cell_type",
    )
    .mark_bar()
)
ch.save(f"{artifact_dir}/nk_and_t_clusters_per_timepoint.svg")
ch.display()

# %% [markdown]
# ## Save results

# %%
adata_nkt.write_h5ad(f"{artifact_dir}/adata_nkt.h5ad")
