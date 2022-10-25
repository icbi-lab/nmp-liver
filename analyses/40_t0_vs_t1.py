# -*- coding: utf-8 -*-
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
adata_nkt_path = nxfvars.get("adata_nkt_path", "../data/results/12_nk_and_t_analysis/artifacts/adata_nkt.h5ad")
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
adata_nkt = sc.read_h5ad(adata_nkt_path)

# %%
de_res = {}
for f in Path(de_res_dir).glob("**/*DESeq2_result.tsv"):
    ct = f.name.replace("_DESeq2_result.tsv", "").split("_", maxsplit=3)[-1]
    de_res[ct] = pd.read_csv(f, sep="\t")

# %%
de_res_nk_and_t = de_res.pop("nk_and_t")

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
stat_mat = (
    pd.concat(
        [
            r.loc[:, ["gene_id", "stat"]]
            .rename(columns={"stat": ct})
            .set_index("gene_id")
            for ct, r in de_res.items()
        ],
        axis=1,
        sort=True,
    )
    .fillna(0)
    .T
)
p_mat.name = "deseq2_stat"

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
pb_nkt = dc.get_pseudobulk(
    adata_nkt,
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
with plt.rc_context({"figure.figsize": (6, 6)}):
    for col in ["cell_type", "cell_type_coarse", "timepoint"]:
        sh.colors.set_scale_anndata(adata_t0t1, col)
        fig = sc.pl.umap(adata_t0t1, color=col, return_fig=True, size=10)
        fig.savefig(f"{artifact_dir}/umap_only_t0t1_{col}.pdf", bbox_inches="tight")

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
cell_type_coarse_fractions = (
    adata_t0t1.obs.groupby(["patient_id", "timepoint"])
    .apply(
        lambda x: x["cell_type_coarse"]
        .value_counts(normalize=True)
        .rename_axis("cell_type_coarse")
    )
    .reset_index(name="frac")
)
tmp_df = (
    cell_type_coarse_fractions.groupby(["timepoint", "cell_type_coarse"])
    .agg("mean")
    .reset_index()
)
ch = (
    alt.Chart(tmp_df)
    .encode(
        x="timepoint",
        y="frac",
        color=alt.Color(
            "cell_type_coarse", scale=sh.colors.altair_scale("cell_type_coarse")
        ),
    )
    .mark_bar()
)
ch.save(
    f"{artifact_dir}/cell_type_coarse_fractions_patient_average_stacked_bar_chart.svg"
)
ch.display()

# %%
ad_ct = sc.AnnData(
    cell_type_fractions.pivot(
        index=["patient_id", "timepoint"], columns="cell_type", values="frac"
    )
)
ad_ct.obs = ad_ct.obs.reset_index()

# %%
fig = sh.pairwise.plot_paired(
    ad_ct,
    "timepoint",
    paired_by="patient_id",
    n_cols=8,
    panel_size=(2.4, 3.5),
    ylabel="frac",
    palette=sh.colors.COLORS.patient_id,
    return_fig=True,
)
fig.savefig(
    f"{artifact_dir}/cell_type_fractions_patient_boxplot.pdf", bbox_inches="tight"
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
ch = sh.compare_groups.pl.plot_lm_result_altair(
    tmp_de_res,
    x="gene_id",
    y="cell_type",
    p_col="padj",
    color="log2FoldChange",
    p_cutoff=np.inf,
)
ch.save(f"{artifact_dir}/t0_vs_v1_interleukins_heatmap.svg")
ch.display()

# %%
fig = sc.pl.dotplot(
    adata_t0t1[adata_t0t1.obs["cell_type"] == "Neutrophils"],
    var_names=gene_set_il["gene_symbol"],
    groupby="timepoint",
    cmap="coolwarm",
    swap_axes=True,
    title="Neutrophils",
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/t0_vs_t1_dotplot_neutro.pdf", bbox_inches="tight")

# %%
fig = sc.pl.dotplot(
    adata_t0t1[adata_t0t1.obs["cell_type"] == "Monocytes ⁄ Macrophages"],
    var_names=gene_set_il["gene_symbol"],
    groupby="timepoint",
    cmap="coolwarm",
    swap_axes=True,
    title="Macro/Mono",
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/t0_vs_t1_dotplot_macro_mono.pdf", bbox_inches="tight")

# %%
fig = sc.pl.dotplot(
    adata_t0t1[adata_t0t1.obs["cell_type_coarse"].isin(["NK cells", "T cells"]), :],
    var_names=gene_set_il["gene_symbol"],
    groupby="timepoint",
    cmap="coolwarm",
    swap_axes=True,
    title="T/NK",
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/t0_vs_t1_dotplot_t_nk.pdf", bbox_inches="tight")

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
ch = sh.compare_groups.pl.plot_lm_result_altair(
    tf_res, y="cell_type", color="act_score", x="TF", title="TFs scores by cell-type"
)
ch.save(f"{artifact_dir}/t0_vs_t1_dorothea_tf_heatmap.svg")
ch.display()

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
msigdb

# %% [markdown]
# ### ORA

# %%
de_genes_per_cell_type = (
    de_res_all.loc[lambda x: (x["padj"] < 0.01) & (np.abs(x["log2FoldChange"] > 1))]
    .groupby("cell_type")
    .size()
    .reset_index(name="n")
)
de_genes_per_cell_type

# %%
median_de_genes = int(de_genes_per_cell_type.agg(("median"))["n"])

# %%
de_res_all.groupby(["cell_type"]).apply(
    lambda x: x.sort_values("pvalue")
    .loc[lambda x: np.abs(x["log2FoldChange"] > 1)]
    .iloc[:median_de_genes]
).reset_index(drop=True).groupby("cell_type").size().reset_index()

# %%
ora_pvals = dc.get_ora_df(
    de_res_all.groupby(["cell_type"])
    .apply(
        lambda x: x.sort_values("pvalue")
        .loc[lambda x: np.abs(x["log2FoldChange"] > 1)]
        .iloc[:median_de_genes]
    )
    .reset_index(drop=True),
    # de_res_all.loc[lambda x: (x["padj"] < 0.01) & (np.abs(x["log2FoldChange"] > 1))],
    msigdb,
    groupby="cell_type",
    features="gene_id",
    source="geneset",
    target="genesymbol",
)
# add pseudocount of 1e-10 to avoid inf value. This means the max "score" is 10
ora_score = -np.log10(ora_pvals + 1e-10)
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
    ),
    how="outer",
).fillna({"act_score": 0, "pvalue": 1, "fdr": 1})

# %%
ch = sh.compare_groups.pl.plot_lm_result_altair(
    ora_res.rename(columns={"act_score": "-log10(p)"}),
    y="cell_type",
    color="-log10(p)",
    x="ORA",
    title="ORA scores by cell-type",
    cmap="viridis",
    domain=lambda x: [0, x],
    reverse=False,
    value_max=10,
    p_cutoff=np.inf,
    cluster=True,
)
ch.save(f"{artifact_dir}/t0_vs_t1_hallmark_ora_fixed_gene_set_size.svg")
ch.display()

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
cpdb_res_n.to_csv(f"{artifact_dir}/cell2cell_res_neutro.csv")

# %%
ch = cpdba.plot_result(
    cpdb_res_n,
    group_col="comparison",
    aggregate=False,
    cluster="heatmap",
    title="CellChat analysis: Neutrophils T0 vs T1",
)
ch.save(f"{artifact_dir}/t0_vs_t1_cellchat_neutro.svg")
ch.display()

# %%
cpdb_res_m = cpdba.significant_interactions(
    de_res["Monocytes ⁄ Macrophages"], max_pvalue=0.1, pvalue_col="pvalue"
)

# %%
cpdb_res_m.to_csv(f"{artifact_dir}/cell2cell_res_myeloid.csv")

# %%
cpdb_res_m

# %%
ch = cpdba.plot_result(
    cpdb_res_m.loc[lambda x: x["log2FoldChange"] > 0],
    group_col="comparison",
    aggregate=False,
    cluster="heatmap",
    title="CellChat analysis: monocytic lineage T0 vs T1",
)
ch.save(f"{artifact_dir}/t0_vs_t1_cellchat_macro_mono.svg")
ch.display()

# %% [markdown]
# ## Neutrophils

# %%
sc.pl.umap(adata_n, color=["patient_id", "timepoint", "cell_type"])

# %%
fig = sc.pl.umap(
    adata_n, color=["CXCR2", "CXCR1", "CXCR4"], cmap="inferno", size=15, return_fig=True
)
fig.savefig(f"{artifact_dir}/neutro_umap_cxcr.pdf", bbox_inches="tight", dpi=600)

# %%
genes = ["CXCR2", "CXCR1", "CXCL8", "CXCR4"]
pvalues = de_res["Neutrophils"].set_index("gene_id").loc[genes, "padj"].tolist()
fig = sh.pairwise.plot_paired(
    pb_n,
    groupby="timepoint",
    paired_by="patient_id",
    var_names=genes,
    pvalues=pvalues,
    pvalue_template=lambda x: f"FDR={x:.4f}, DESeq2"
    if x >= 0.0001
    else "FDR<0.0001, DESeq2",
    palette=sh.colors.COLORS.patient_id,
    panel_size=(2.5, 4),
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/neutro_t0_vs_t1_boxplot_cxcr.pdf", bbox_inches="tight")

# %%
for gene in genes:
    fig = sc.pl.umap(adata, color=gene, cmap="inferno", return_fig=True)
    fig.savefig(f"{artifact_dir}/umap_{gene}.pdf", dpi=600)

# %%
neutro_goi = "CXCR1, CXCR2, CXCR4, CXCL8, VEGFA, ARG1, CD83, ICAM1, PTGS2, SELL, CCRL2, CCL3, CCL4".replace(
    " ", ""
).split(
    ","
)

# %%
for name, genes in {"neutro_genes_of_interest": neutro_goi}.items():
    pvalues = (
        de_res["Neutrophils"]
        .set_index("gene_id")
        .loc[genes, "padj"]
        .tolist()
    )
    fig = sh.pairwise.plot_paired(
        pb_n,
        groupby="timepoint",
        paired_by="patient_id",
        var_names=genes,
        pvalues=pvalues,
        pvalue_template=lambda x: f"{x:.3f}" if x >= 0.001 else "<0.001",
        palette=sh.colors.COLORS.patient_id,
        return_fig=True,
        size=8,
        panel_size=(1, 4),
        n_cols=15,
        show=False,
    )
    for i, ax in enumerate(fig.axes):
        ax.set_ylim(
            np.min(pb_n[:, genes].X) - 0.5,
            np.max(pb_n[:, genes].X) + 0.5,
        )
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        if i > 0:
            ax.yaxis.set_ticklabels([])
            ax.set_ylabel(None)
    plt.subplots_adjust(wspace=0.25)
    fig.savefig(f"{artifact_dir}/neutro_t0_vs_t1_boxplot_tight_{name}.pdf", bbox_inches="tight")

# %%
ch = sh.pairwise.plot_paired_fc(
    pb_n,
    groupby="timepoint",
    paired_by="patient_id",
    var_names=neutro_goi,
).properties(height=200)
ch.save(f"{artifact_dir}/neutro_t0_vs_t1_fold_change_barchart_neutro_genes_of_interest.svg")
ch.display()

# %%
top_genes = de_res["Neutrophils"].assign(direction = lambda x: np.sign(x["log2FoldChange"])).groupby("direction").apply(lambda x: x.iloc[:30])["gene_id"]
ch = sh.pairwise.plot_paired_fc(
    pb_n,
    groupby="timepoint",
    paired_by="patient_id",
    var_names=top_genes,
).properties(height=200)
ch.save(f"{artifact_dir}/neutro_t0_vs_t1_fold_change_barchart_top_genes.svg")
ch.display()

# %%
fig = sc.pl.matrixplot(pb_n, var_names=top_genes, cmap="viridis", groupby="timepoint", return_fig=True)
fig.savefig(f"{artifact_dir}/neutro_t0_vs_t1_heatmap_top_genes.pdf", bbox_inches="tight")

# %% [markdown]
# ## monocytic lineage

# %%
sc.pl.umap(adata_m, color=["patient_id", "timepoint", "cell_type"])

# %%
genes = ["LYZ", "FCN1", "VCAN", "HLA-DRA", "CD163", "MARCO", "HMOX1", "VSIG4"]
pvalues = (
    de_res["Monocytes ⁄ Macrophages"].set_index("gene_id").loc[genes, "padj"].tolist()
)
fig = sh.pairwise.plot_paired(
    pb_m,
    groupby="timepoint",
    paired_by="patient_id",
    var_names=genes,
    pvalues=pvalues,
    pvalue_template=lambda x: f"FDR={x:.4f}, DESeq2"
    if x >= 0.0001
    else "FDR<0.0001, DESeq2",
    palette=sh.colors.COLORS.patient_id,
    return_fig=True,
    panel_size=(2.5, 4),
)
fig.savefig(f"{artifact_dir}/myeloid_t0_vs_t1_boxplot.pdf", bbox_inches="tight")

# %%
goi_m_inflammatory = "LYZ, FCN1, VCAN, HLA-DRA, S100A8, S100A9, S100A12, CXCL8, MNDA, CSTA, CD74, CCL2, CCL3".replace(
    " ", ""
).split(
    ","
)
goi_m_tolerogenic = (
    " CD163, MARCO, HMOX1, VSIG4, CD5L, NSMAF, CTSB, VMO1, HMOX1".replace(
        " ", ""
    ).split(",")
)

# %%
for name, genes in {"inflammatory": goi_m_inflammatory, "tolerogenic": goi_m_tolerogenic}.items():
    pvalues = (
        de_res["Monocytes ⁄ Macrophages"]
        .set_index("gene_id")
        .loc[genes, "padj"]
        .tolist()
    )
    fig = sh.pairwise.plot_paired(
        pb_m,
        groupby="timepoint",
        paired_by="patient_id",
        var_names=genes,
        pvalues=pvalues,
        pvalue_template=lambda x: f"{x:.3f}" if x >= 0.001 else "<0.001",
        palette=sh.colors.COLORS.patient_id,
        return_fig=True,
        size=8,
        panel_size=(1, 4),
        n_cols=15,
        show=False,
    )
    for i, ax in enumerate(fig.axes):
        ax.set_ylim(
            np.min(pb_m[:, genes].X) - 0.5,
            np.max(pb_m[:, genes].X) + 0.5,
        )
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
        if i > 0:
            ax.yaxis.set_ticklabels([])
            ax.set_ylabel(None)
    plt.subplots_adjust(wspace=0.25)
    fig.savefig(f"{artifact_dir}/myeloid_t0_vs_t1_boxplot_tight_{name}.pdf", bbox_inches="tight")

# %%
for name, genes in {"inflammatory": goi_m_inflammatory, "tolerogenic": goi_m_tolerogenic}.items():
    ch = sh.pairwise.plot_paired_fc(
        pb_m,
        groupby="timepoint",
        paired_by="patient_id",
        var_names=genes,
    ).properties(height=200)
    ch.save(f"{artifact_dir}/myeloid_t0_vs_t1_fold_change_barchart_{name}.svg")
    ch.display()

# %%
top_genes = de_res["Monocytes ⁄ Macrophages"].assign(direction = lambda x: np.sign(x["log2FoldChange"])).groupby("direction").apply(lambda x: x.iloc[:30])["gene_id"]
ch = sh.pairwise.plot_paired_fc(
    pb_m,
    groupby="timepoint",
    paired_by="patient_id",
    var_names=top_genes,
).properties(height=200)
ch.save(f"{artifact_dir}/myeloid_t0_vs_t1_fold_change_barchart_top_genes.svg")
ch.display()

# %%
fig = sc.pl.matrixplot(pb_m, var_names=top_genes, cmap="viridis", groupby="timepoint", return_fig=True)
fig.savefig(f"{artifact_dir}/myeloid_t0_vs_t1_heatmap_top_genes.pdf", bbox_inches="tight")

# %% [markdown]
# ## NK and T cells

# %%
top_genes = de_res_nk_and_t.assign(direction = lambda x: np.sign(x["log2FoldChange"])).groupby("direction").apply(lambda x: x.iloc[:30])["gene_id"]
ch = sh.pairwise.plot_paired_fc(
    pb_nkt,
    groupby="timepoint",
    paired_by="patient_id",
    var_names=top_genes,
).properties(height=200)
ch.save(f"{artifact_dir}/nk_and_t_t0_vs_t1_fold_change_barchart_top_genes.svg")
ch.display()

# %% [markdown]
# ## QC plots

# %%
df_qc = adata_t0t1.obs.groupby(["patient_id", "timepoint"], observed=True).agg(
    umi_counts=("total_counts", "mean"),
    n_genes=("n_genes_by_counts", "mean"),
    pct_mito=("pct_counts_mito", "mean"),
).reset_index()
tmp_ad = sc.AnnData(df_qc.loc[:, ["umi_counts", "n_genes", "pct_mito"]], obs=df_qc)

# %%
tmp_ad

# %%
fig = sh.pairwise.plot_paired(tmp_ad, "timepoint", paired_by="patient_id", palette=sh.colors.COLORS.patient_id, return_fig=True)
fig.savefig(f"{artifact_dir}/t0_vs_t1_qc_metrics.pdf", bbox_inches="tight")

# %%
qc_genes=["FOS", "JUN", "JUNB", "JUND", "ATF3", "EGR1", "HSPA1A", "HSPA1B", "HSP90AB1", "HSPA8", "IER3", "IER2", "BTG1", "BTG2", "DUSP1", "LMNA", "CYCS", "DDIT3", "APAF1", "BAX", "BAK1", "CASP3", "CASP9"]

# %%
fig = sc.pl.dotplot(adata_t0t1, var_names=qc_genes, groupby="timepoint", cmap="viridis", return_fig=True)
fig.savefig(f"{artifact_dir}/t0_vs_t1_stress_apoptosis_genes_dotplot.pdf", bbox_inches="tight")

# %%
