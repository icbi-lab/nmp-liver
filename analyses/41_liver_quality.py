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
from itertools import zip_longest
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
    adata_t0t1,
    groupby=["patient_id", "timepoint", "cell_type", "LT", "ECD"],
    log_norm=True,
)

# %% [markdown]
# ## Comparison by quality and timepoints

# %%
df_neutro = (
    pb[:, gene_set_il["gene_symbol"]]
    .to_df()
    .join(pb.obs)
    .loc[lambda x: x["cell_type"] == "Neutrophils"]
)

# %%
df_m = (
    pb[:, gene_set_il["gene_symbol"]]
    .to_df()
    .join(pb.obs)
    .loc[lambda x: x["cell_type"] == "monocytic lineage"]
)


# %%
def make_paired_plot(tmp_df, hue):
    var_names = gene_set_il["gene_symbol"]
    fig, axes = plt.subplots(3, 10, figsize=(20, 10), squeeze=False)
    axes = axes.flatten()
    max_val = np.max(tmp_df.loc[:, var_names].values)
    for i, (var, ax) in enumerate(zip_longest(var_names, axes)):
        if var is not None:
            sns.stripplot(
                x="timepoint",
                data=tmp_df,
                y=var,
                ax=ax,
                hue=hue,
                size=5,
                linewidth=1,
                palette=getattr(sh.colors.COLORS, hue),
                dodge=True,
            )
            sns.boxplot(
                tmp_df, x="timepoint", hue=hue, y=var, ax=ax, color="white", fliersize=0
            )
            sns.lineplot(
                tmp_df,
                x="timepoint",
                hue=hue,
                y=var,
                errorbar=None,
                legend=False,
                palette=getattr(sh.colors.COLORS, hue),
                linewidth=3,
                ax=ax,
            )
            ax.legend().set_visible(False)
            ax.set_ylim(-0.5, max_val)
            ax.set_title(var)

        else:
            ax.set_visible(False)

    fig.tight_layout()
    sh.util.adjust_box_widths(fig, 0.8)
    # return fig


# %%
make_paired_plot(df_neutro, "LT")

# %%
make_paired_plot(df_neutro, "ECD")

# %%
make_paired_plot(df_m, "LT")

# %%
make_paired_plot(df_m, "ECD")

# %%
