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
adata_neutro_path = nxfvars.get(
    "adata_neutro_path", "../data/results/11_neutro_analysis/artifacts/adata_n.h5ad"
)
adata_myeloid_path = nxfvars.get(
    "adata_myeloid_path", "../data/results/10_myeloid_analysis//artifacts/adata_m.h5ad"
)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/nmp-temp/")
de_res_dir = nxfvars.get("de_res_dir", "../data/results/21_deseq/DESEQ_T0_T1/")
dorothea_net = nxfvars.get("dorothea", "../tables/dorothea_human_AB_2022-09-28.csv")
cellchatdb_path = nxfvars.get("cellchatdb", "../tables/cellchatdb_2022-09-29.tsv")
msigdb_path = nxfvars.get("msigdb", "../tables/gene_sets_hallmarks_msigdb.csv")
gene_set_il_path = nxfvars.get(
    "gene_set_il_path", "../tables/gene_sets_interleukins_chemokines.xlsx"
)

# %% [markdown]
# ## Load data

# %%
adata = sc.read_h5ad(adata_path)

# %%
tfnet = pd.read_csv(dorothea_net)

# %%
cellchatdb = pd.read_csv(cellchatdb_path, sep="\t", comment="#")

# %%
msigdb = pd.read_csv(msigdb_path, comment="#").drop_duplicates()

# %%
gene_set_il = pd.read_excel(gene_set_il_path)

# %%
adata_t0t1 = adata[adata.obs["timepoint"].isin(["T0", "T1"]), :]

# %%
adata_n = sc.read_h5ad(adata_neutro_path)
adata_m = sc.read_h5ad(adata_myeloid_path)

# %%
de_res = {}
for f in Path(de_res_dir).glob("**/*.tsv"):
    ct = f.name.replace("_DESeq2_result.tsv", "").split("_")[-1]
    de_res[ct] = pd.read_csv(f, sep="\t")

# %%
de_res_all = pd.concat([df.assign(cell_type=ct) for ct, df in de_res.items()])

# %%
logfc_mat = (
    pd.concat(
        [
            r.loc[:, ["gene_id", "log2FoldChange"]]
            .rename(columns={"log2FoldChange": ct})
            .set_index("gene_id")
            for ct, r in de_res.items()
        ],
        axis=1,
        sort=True,
    )
    .fillna(0)
    .T
)
logfc_mat.name = "deseq2_estimate"

# %%
p_mat = (
    pd.concat(
        [
            r.loc[:, ["gene_id", "padj"]]
            .rename(columns={"padj": ct})
            .set_index("gene_id")
            for ct, r in de_res.items()
        ],
        axis=1,
        sort=True,
    )
    .fillna(1)
    .T
)
p_mat.name = "deseq2_pvals"

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
sc.pp.normalize_total(pb, target_sum=1e6)
sc.pp.log1p(pb)

# %%
pb_n = dc.get_pseudobulk(
    adata_n,
    sample_col="patient_id",
    groups_col="timepoint",
    min_prop=0.05,
    min_cells=10,
    min_counts=1000,
    min_smpls=2,
)
sc.pp.normalize_total(pb_n, target_sum=1e6)
sc.pp.log1p(pb_n)

# %%
pb_m = dc.get_pseudobulk(
    adata_m,
    sample_col="patient_id",
    groups_col="timepoint",
    min_prop=0.05,
    min_cells=10,
    min_counts=1000,
    min_smpls=2,
)
sc.pp.normalize_total(pb_m, target_sum=1e6)
sc.pp.log1p(pb_m)

# %%
sc.pl.umap(adata, color="cell_type")

# %% [markdown]
# ## Overview

# %%
sc.pl.umap(adata_t0t1, color="cell_type")

# %%
sh.colors.set_scale_anndata(adata_t0t1, "timepoint")

# %%
sc.pl.umap(adata_t0t1, color="timepoint", size=15)

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
# ### Highlighted genes

# %%
# Fill genes with logFC=0, FDR=1 if they are not contained for a certain cell-type
tmp_de_res = (
    de_res_all.merge(
        pd.DataFrame.from_records(
            itertools.product(
                de_res_all["gene_id"].unique(), de_res_all["cell_type"].unique()
            ),
            columns=["gene_id", "cell_type"],
        ),
        how="right",
        on=["gene_id", "cell_type"],
    )
    .assign(
        log2FoldChange=lambda x: x["log2FoldChange"].fillna(0),
        padj=lambda x: x["padj"].fillna(1),
    )
    .loc[lambda x: x["gene_id"].isin(gene_set_il["gene_symbol"])]
)

# %%
sh.compare_groups.pl.plot_lm_result_altair(
    tmp_de_res,
    x="gene_id",
    y="cell_type",
    p_col="padj",
    color="log2FoldChange",
    p_cutoff=np.inf
)

# %%
sc.pl.dotplot(
    adata_t0t1[adata_t0t1.obs["cell_type"] == "Neutrophils"],
    var_names=gene_set_il["gene_symbol"],
    groupby="timepoint",
    cmap="coolwarm",
    swap_axes=True,
    title="Neutrophils",
)

# %%
sc.pl.dotplot(
    adata_t0t1[adata_t0t1.obs["cell_type"] == "monocytic lineage"],
    var_names=gene_set_il["gene_symbol"],
    groupby="timepoint",
    cmap="coolwarm",
    swap_axes=True,
    title="Macro/Mono",
)

# %% [markdown]
# ## Transcription factors (Dorothea)

# %%
# TFs model was obtained with this function from omnipathdb.
# For reproducibility, the result was stored in the `tables` directory and is used for all analyses.
# tfnet = dc.get_dorothea(organism="human", levels=["A", "B"])

# %%
tf_acts, tf_pvals = dc.dense_run(dc.run_mlm, mat=logfc_mat, net=tfnet, verbose=True)


# %%
def format_decoupler_results(tf_acts, tf_pvals, name="variable", contrast="contrast"):
    return (
        dc.format_contrast_results(tf_acts, tf_pvals)
        .drop(columns="adj_pvals")
        .rename(
            columns={
                "logFCs": "act_score",
                "pvals": "pvalue",
                "name": name,
                "contrast": contrast,
            }
        )
        .assign(
            act_score=lambda x: x["act_score"].fillna(0),
            pvalue=lambda x: x["pvalue"].fillna(1),
        )
        .pipe(sh.util.fdr_correction)
    )


# %%
tf_res = format_decoupler_results(tf_acts, tf_pvals, name="TF", contrast="cell_type")

# %%
sh.compare_groups.pl.plot_lm_result_altair(
    tf_res, y="cell_type", color="act_score", x="TF", title="TFs scores by cell-type"
)

# %% [markdown]
# ### Check individual pathway 
#  * using PPARG as an example
#  * check if the results make sense 
#  * check if the direction of the score makes sense (red = up) 

# %%
dc.plot_volcano(logfc_mat, p_mat, "Neutrophils", name="PPARG", net=tfnet, top=5)

# %%
sh.pairwise.plot_paired(
    pb,
    "timepoint",
    var_names=["ACSL5", "CDKN1A", "CD83", "NR1D1", "PLIN2"],
    paired_by="patient_id",
)

# %% [markdown]
# ## GSEA

# %%
# # MSigDB was obtained as follows. The results were stored in `tables` for reproducibility reasons
# msigdb = dc.get_resource('MSigDB')
# genesets = pd.read_csv("../tables/gene_sets_hallmarks.csv", header=None)[0].tolist()
# for gs in genesets:
#     assert np.sum(msigdb["geneset"] == gs), gs
# msigdb_subset = msigdb.loc[lambda x: (x["collection"] == "hallmark") & (x["geneset"].isin(genesets)), :]
# msigdb_subset.to_csv("../tables/gene_sets_hallmarks_msigdb.csv", index=False)

# %%

# %%

# %%
msigdb

# %%
ora_pvals = dc.get_ora_df(
    de_res_all.loc[lambda x: (x["padj"] < 0.01) & (np.abs(x["log2FoldChange"] > 1))],
    msigdb,
    groupby="cell_type",
    features="gene_id",
    source="geneset",
    target="genesymbol",
)
ora_score = -np.log10(ora_pvals)
ora_score.name = "ora_estimate"

# %%
ora_res = format_decoupler_results(
    ora_score, ora_pvals, name="ORA", contrast="cell_type"
)

# %%
# re-add cell-types that were removed because they have no significant genes
ora_res = ora_res.merge(
    pd.DataFrame.from_records(
        itertools.product(de_res_all["cell_type"].unique(), ora_res["ORA"].unique()),
        columns=["cell_type", "ORA"],
    ), how="outer"
).fillna({"act_score": 0, "pvalue": 1, "fdr": 1})

# %%
sh.compare_groups.pl.plot_lm_result_altair(
    ora_res.rename(columns={"act_score": "-log10(p)"}),
    y="cell_type",
    color="-log10(p)",
    x="ORA",
    title="TFs scores by cell-type",
    cmap="viridis",
    domain=lambda x: [0, x],
    reverse=False,
    value_max=10,
)

# %% [markdown]
# ## Cellphonedb analysis

# %%
cpdba = sh.cell2cell.CpdbAnalysis(
    cellchatdb,
    adata_t0t1,
    pseudobulk_group_by=["patient_id"],
    cell_type_column="cell_type",
)

# %%
cpdb_res_n = cpdba.significant_interactions(
    de_res["Neutrophils"], max_pvalue=0.1, pvalue_col="pvalue"
)

# %%
cpdba.plot_result(
    cpdb_res_n,
    group_col="comparison",
    aggregate=False,
    cluster="heatmap",
    title="CellChat analysis: Neutrophils T0 vs T1",
)

# %%
cpdb_res_m = cpdba.significant_interactions(
    de_res["monocytic lineage"], max_pvalue=0.1, pvalue_col="pvalue"
)

# %%
cpdba.plot_result(
    cpdb_res_m,
    group_col="comparison",
    aggregate=False,
    cluster="heatmap",
    title="CellChat analysis: monocytic lineage T0 vs T1",
)

# %% [markdown]
# ## Neutrophils

# %%
sc.pl.dotplot(
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
)

# %%
sc.pl.umap(adata_n, color=["patient_id", "timepoint", "cell_type"])

# %%
sc.pl.umap(adata_n, color=["CXCR2", "CXCR1", "CXCR4"], cmap="inferno", size=15)

# %%
genes = ["CXCR2", "CXCR1", "CXCL8", "CXCR4"]
pvalues = de_res["Neutrophils"].set_index("gene_id").loc[genes, "padj"].tolist()
sh.pairwise.plot_paired(
    pb_n[pb_n.obs["timepoint"].isin(["T0", "T1"]), :],
    groupby="timepoint",
    paired_by="patient_id",
    var_names=genes,
    pvalues=pvalues,
    pvalue_template=lambda x: f"FDR={x:.4f}, DESeq2"
    if x >= 0.0001
    else "FDR<0.0001, DESeq2",
    palette=sh.colors.COLORS.patient_id,
)

# %%
sc.pl.umap(adata, color="CXCL8", cmap="inferno")

# %% [markdown]
# ## monocytic lineage

# %%
sc.pl.umap(adata_m, color=["patient_id", "timepoint", "cell_type"])

# %%
sc.pl.dotplot(
    adata,
    groupby="cell_type_coarse",
    var_names=[
        "CST3",
        "CSTB",
        "MS4A7",
        "MARCH1",
        "HMOX1",
        "CD68",
        "MAFB",
        "CD163",
        "VCAN",
        "CSF1R",
        "CD300E",
        "LIPA",
        "SAMHD1",
        "PSAP",
        "TGFBI",
    ],
    cmap="coolwarm",
)

# %%
genes = ["LYZ", "FCN1", "VCAN", "HLA-DRA", "CD163", "MARCO", "HMOX1", "VSIG4"]
pvalues = de_res["monocytic lineage"].set_index("gene_id").loc[genes, "padj"].tolist()
sh.pairwise.plot_paired(
    pb_m[pb_m.obs["timepoint"].isin(["T0", "T1"]), :],
    groupby="timepoint",
    paired_by="patient_id",
    var_names=genes,
    pvalues=pvalues,
    pvalue_template=lambda x: f"FDR={x:.4f}, DESeq2"
    if x >= 0.0001
    else "FDR<0.0001, DESeq2",
    palette=sh.colors.COLORS.patient_id,
)

# %%
