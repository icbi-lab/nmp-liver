# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: 'SSH apollo-15 apollo-15: 2022-schneeberger-liver'
#     language: ''
#     name: rik_ssh_apollo_15_apollo152022schneebergerliver
# ---

# %%
import scvi
import scanpy as sc
import pandas as pd
from nxfvars import nxfvars
from threadpoolctl import threadpool_limits
from tqdm.contrib.concurrent import process_map
import multiprocessing

# %%
sc.settings.set_figure_params(figsize=(4,4), frameon=False)


# %%
def set_all_seeds(seed=0):
    import os
    import random
    import numpy as np
    import torch

    scvi.settings.seed = seed
    os.environ["PYTHONHASHSEED"] = str(seed)  # Python general
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    np.random.seed(seed)  # Numpy random
    random.seed(seed)  # Python random

    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # For multiGPU


# %%
adata_path = nxfvars.get("adata_path", "../data/results/02_qc_and_filtering/NMP_Liver.qc.h5ad")
cpus = nxfvars.get("cpus", 8)
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads/nmp-temp")

# %%
threadpool_limits(cpus)

# %%
set_all_seeds()

# %%
adata = sc.read_h5ad(adata_path)

# %% [markdown]
# ## Highly variable genes
#
# "Between 1000 and 10000" standard for scvi analysis. 
# For the lung atlas I used 6000 and it was a bit too much (sometimes imperfect integration), so let's try a bit less this time. 

# %%
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=4000,
    subset=False,
    flavor="seurat_v3",
    batch_key="patient_id"
)

# %%
pd.set_option("display.max_columns", None)

# %% [markdown]
# ## Run SCVI

# %%
adata_scvi = adata[:, adata.var["highly_variable"]].copy()

# %%
scvi.model.SCVI.setup_anndata(adata_scvi, batch_key="sample_id")

# %%
model = scvi.model.SCVI(adata_scvi)

# %%
model.train(early_stopping=True)

# %% [markdown]
# ## Standard scanpy preprocessing (for comparison)

# %%
adata_raw = adata.copy()
sc.pp.normalize_total(adata_raw)
sc.pp.log1p(adata_raw)
sc.tl.pca(adata_raw, use_highly_variable=True)
adata.raw =adata_raw

# %%
adata.obsm["X_pca"] = adata_raw.obsm["X_pca"]

# %%
sc.pp.neighbors(adata, key_added="neighbors_uncorrected")

# %%
sc.tl.umap(adata, neighbors_key="neighbors_uncorrected")
adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"]
del adata.obsm["X_umap"]

# %%
adata.obs

# %%
sc.pl.embedding(adata, "umap_uncorrected", color=["patient_id", "timepoint", "FCGR3B", "CD68", "CD3E", "ALB", "n_genes", "pct_counts_mito"], cmap="inferno")

# %% [markdown]
# ## With batch correction

# %%
adata.obsm["X_scVI"] = model.get_latent_representation()

# %%
sc.pp.neighbors(adata, use_rep="X_scVI")

# %%
sc.tl.umap(adata)

# %%
sc.pl.umap(adata, color=["patient_id", "timepoint", "n_genes", "pct_counts_mito"])

# %%
sc.pl.umap(adata, color=["FCGR3B", "CD68", "FLT3", "CD3E", "GNLY", "FLT1", "CEACAM8", "ALB", "KRT7", "CD79A", "MZB1", "B2M"], cmap="inferno")

# %% [markdown]
# ## Doublet detection (SOLO)

# %%
solo_models = [scvi.external.SOLO.from_scvi_model(model, restrict_to_batch=batch) for batch in adata.obs["sample_id"].unique()]


# %%
def run_solo(solo_model):
    solo_model.train()
    return solo_model.predict(soft=False)


# %%
solo_res = []
for solo_model in solo_models:
    solo_res.append(run_solo(solo_model))

# %%
adata.obs["is_doublet"] = pd.concat(solo_res)

# %%
solo_series = pd.concat(solo_res)
solo_series.index = solo_series.index.str.replace("-0$", "", regex=True)
adata.obs["is_doublet"] = solo_series

# %%
sc.pl.umap(adata, color="is_doublet")

# %% [markdown]
# ## Filter doublets and reprocess

# %%
adata_nodoublet = adata[adata.obs["is_doublet"] == "singlet", :].copy()
adata_nodoublet.obsm["X_scVI"] = model.get_latent_representation(adata_scvi[adata_nodoublet.obs_names, :])

# %%
sc.pp.neighbors(adata_nodoublet, use_rep="X_scVI")

# %%
sc.tl.umap(adata_nodoublet)

# %%
sc.pl.umap(adata_nodoublet, color=["patient_id", "timepoint", "n_genes", "pct_counts_mito"])

# %%
sc.pl.umap(adata_nodoublet, color=["FCGR3B", "CD68", "FLT3", "CD3E", "GNLY", "FLT1", "CEACAM8", "ALB", "KRT7", "CD79A", "MZB1", "B2M"], cmap="inferno")

# %% [markdown]
# ## Save results

# %%
model.save(f"{artifact_dir}/scvi_model", save_anndata=True)

# %%
adata.write_h5ad(f"{artifact_dir}/adata_scvi_before_doublets.h5ad")

# %%
adata_nodoublet.write_h5ad(f"{artifact_dir}/adata_scvi_doublet_filtered.h5ad")

# %%
adata.shape

# %%
adata_nodoublet.shape

# %%
