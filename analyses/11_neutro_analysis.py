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
import scanpy as sc
import scanpy_helpers as sh
import pandas as pd
from nxfvars import nxfvars
import matplotlib.pyplot as plt
import altair as alt
from tqdm import tqdm

# %%
sc.settings.set_figure_params(figsize=(4, 4), frameon=False)

# %%
adata_path = nxfvars.get(
    "adata_path", "../data/results/04_cell_type_annotation/artifacts/adata_cell_types.h5ad"
)
scvi_model = nxfvars.get("scvi_model", "/home/sturm/Downloads/nmp-temp/scvi_model/")
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/nmp-temp/")

# %%
# model = scvi.model.SCVI.load(scvi_model)

# %%
adata = sc.read_h5ad(adata_path)

# %%
sc.pl.umap(adata, color="cell_type")

# %% [markdown]
# ## Cluster overview

# %%
adata_n = adata[adata.obs["cell_type"].str.startswith("Neutro"), :].copy()

# %%
pb_n = sh.pseudobulk.pseudobulk(adata_n, groupby=["patient_id", "cell_type"])
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n, base=2)

# %%
sc.pp.neighbors(adata_n, use_rep="X_scVI")

# %%
sc.tl.umap(adata_n)

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    sc.pl.umap(adata_n, color="cell_type", legend_loc="on data", legend_fontoutline=2)

# %%
sc.pl.umap(adata_n, color=["patient_id", "timepoint"])

# %%
tmp_df = (
    adata_n.obs.groupby(["patient_id", "cell_type"]).size().reset_index(name="n_cells")
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
(heatmp + txt).properties(width=200)

# %%
tmp_df = (
    adata_n.obs.groupby(["patient_id", "cell_type", "timepoint"])
    .size()
    .reset_index(name="n_cells")
    .groupby(["cell_type", "timepoint"])
    .agg("mean")
    .reset_index()
)

# %%
alt.Chart(tmp_df).encode(
    color="cell_type", y="timepoint", x=alt.X("n_cells", stack="normalize")
).mark_bar()

# %%
alt.Chart(tmp_df).encode(
    color="timepoint", y=alt.Y("n_cells", stack="normalize"), x="cell_type"
).mark_bar()

# %% [markdown]
# ## Characterization of clusters

# %%
for ct in tqdm(pb_n.obs["cell_type"].unique()):
    sh.signatures.fold_change(
        pb_n, obs_col="cell_type", positive_class=ct, key_added=f"{ct}_fc"
    )
    sh.signatures.specific_fold_change(
        pb_n, obs_col="cell_type", positive_class=ct, key_added=f"{ct}_sfc"
    )
    sh.signatures.roc_auc(
        pb_n, obs_col="cell_type", positive_class=ct, key_added=f"{ct}_auroc"
    )

# %%
markers = {
    ct: pb_n.var.loc[
        lambda x: (x[f"{ct}_auroc"] >= 0.7)
        & (x[f"{ct}_fc"] > 2)
        & (x[f"{ct}_sfc"] > 1.5)
    ]
    .sort_values(f"{ct}_auroc", ascending=False)
    .index.tolist()
    for ct in sorted(pb_n.obs["cell_type"].unique())
}

# %%
sh.signatures.plot_markers(pb_n, "cell_type", markers, top=10, return_fig=False)

# %%
sh.signatures.plot_metric_strip(pb_n, markers, top=10)

# %%
with plt.rc_context({"figure.figsize": (6, 6)}):
    sc.pl.umap(adata_n, color="cell_type", legend_loc="on data", legend_fontoutline=2)

# %%
for cluster, m in markers.items():
    print(cluster)
    sc.pl.umap(adata_n, color=m[:10], ncols=5, cmap="inferno", size=10)

# %%