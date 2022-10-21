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
# model = scvi.model.SCVI.load(scvi_model)

# %%
adata = sc.read_h5ad(adata_path)

# %%
sc.pl.umap(adata, color="cell_type")

# %% [markdown]
# ## Markers

# %%
fig = sc.pl.dotplot(
    adata,
    groupby="cell_type_coarse",
    var_names=[
        "IFITM2",
        "CSF3R",
        "H3F3A",
        "S100A11",
        "FPR1",
        "FCGR3B",
        "NAMPT",
        "VNN2",
        "BASP1",
        "G0S2",
        "MXD1",
        "LITAF",
        "CXCR2",
        "S100A9",
        "SOD2",
        "SRGN",
    ],
    cmap="coolwarm",
    return_fig=True
)
fig.savefig(f"{artifact_dir}/neutro_dotplot.pdf", bbox_inches="tight")

# %% [markdown]
# ## Cluster overview

# %%
adata_n = adata[adata.obs["cell_type"].str.startswith("Neutro"), :].copy()

# %%
adata_n.shape

# %%
# Only T0 and T1
adata_n = adata_n[adata_n.obs["timepoint"] != "T2", :].copy()

# %%
adata_n.shape

# %%
sc.pp.neighbors(adata_n, use_rep="X_scVI", n_neighbors=40)

# %%
sc.tl.umap(adata_n)

# %%
sc.tl.leiden(adata_n, resolution=0.25)

# %%
adata_n.obs["cell_type"] = ["N" + x for x in adata_n.obs["leiden"]]

# %%
pb_n = sh.pseudobulk.pseudobulk(adata_n, groupby=["patient_id", "cell_type"])
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n, base=2)

# %%
for col in ["cell_type", "patient_id", "timepoint", "sample_id"]:
    if col != "sample_id":
        sh.colors.set_scale_anndata(
            adata_n, col, "neutro_clusters" if col == "cell_type" else None
        )
    legend_params = (
        dict(legend_loc="on data", legend_fontoutline=2) if col == "cell_type" else {}
    )
    with plt.rc_context({"figure.figsize": (6, 6)}):
        fig = sc.pl.umap(adata_n, color=col, return_fig=True, size=20, **legend_params)
        fig.savefig(
            f"{artifact_dir}/umap_neutrophil_cluster_overview_{col}.pdf",
            bbox_inches="tight", dpi=1200
        )

# %%
adata_n.obs["cell_type"].value_counts()

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
ch = (heatmp + txt).properties(width=200)
ch.save(f"{artifact_dir}/neutrophil_cluster_count_per_patient.svg")
ch.display()

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
ch = (
    alt.Chart(tmp_df)
    .encode(
        color=alt.Color("cell_type", scale=sh.colors.altair_scale("neutro_clusters")),
        y="timepoint",
        x=alt.X("n_cells", stack="normalize"),
    )
    .mark_bar()
)
ch.save(f"{artifact_dir}/timepoints_by_neutro_clusters.svg")
ch.display()

# %%
ch = alt.Chart(tmp_df).encode(
        color=alt.Color("timepoint", scale=sh.colors.altair_scale("timepoint")),
    y=alt.Y("n_cells", stack="normalize"),
    x="cell_type",
).mark_bar()
ch.save(f"{artifact_dir}/neutro_clusters_per_timepoint.svg")

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
        lambda x: (x[f"{ct}_auroc"] >= 0.7) & (x[f"{ct}_fc"] > 1) & (x[f"{ct}_sfc"] > 0)
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
    if len(m):
        sc.pl.umap(adata_n, color=m[:10], ncols=5, cmap="inferno", size=10)

# %% [markdown]
# ## Save results

# %%
with open(f"{artifact_dir}/markers_neutro.csv", "w") as f:
    for ct, genes in markers.items():
        for g in genes:
            f.write(f"{ct},{g}\n")

# %%
adata_n.write_h5ad(f"{artifact_dir}/adata_n.h5ad")
