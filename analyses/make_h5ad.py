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
import pandas as pd
from pathlib import Path
import scanpy as sc
from tqdm.contrib.concurrent import process_map
import scipy.sparse
import anndata

# %% [markdown]
# Abbreviations: LT: liver transplantation, NMP: normothermic machine perfusion, y: years, CVA: cerebrovascular accident, ECD: extended criteria donor, DRI: Donor Risk Index, Indication for NMP (R: recipient, D: donor, L: logistic), CIT: cold ischemia time, PT: preservation time, h: hours
#

# %% [markdown]
# ## Process metadata

# %%
samplesheet = (
    pd.read_excel("../tables/samplesheet_scrnaseq_cleaned.xlsx")
    .rename(
        columns={
            "Liver": "sample_id",
            "Ansatz": "library_id",
            "Datum": "library_date",
            "NGS": "sequencing_id",
        }
    )
    .assign(
        sample_id = lambda x: x["sample_id"].str.replace(" ", "_"),
        library_date=lambda x: pd.to_datetime(x["library_date"], format="%y%m%d").astype(str),
        patient_id=lambda x: x["sample_id"].str.extract(r".*(P\d+).*"),
        timepoint=lambda x: x["sample_id"].str.extract(r".*_(T\d+)$")
    )
)

# %%
n_to_patient = {
    27: "P1",
    28: "P2",
    29: "P3",
    30: "P4",
    31: "P5", 
    32: "P6",
    33: "P7", 
    34: "P8"
}

# %%
patient_meta = pd.read_excel("../tables/patient_metadata.xlsx").assign(patient_id=lambda x: x["n"].map(n_to_patient))

# %%
full_meta = samplesheet.merge(patient_meta, on="patient_id").drop(columns=["n"])

# %%
full_meta

# %%
full_meta.to_excel("/home/sturm/Downloads/samplesheet_full.xlsx")


# %% [markdown]
# ## Load counts

# %%
def load_counts(patient_meta):
    filename = f"../data/01_counts/{patient_meta['patient_id']}{patient_meta['timepoint']}/{patient_meta['sequencing_id']}_RSEC_MolsPerCell.csv"
    gene_expr = pd.read_csv(filename, comment="#", index_col=0)
    adata = sc.AnnData(X=gene_expr)
    adata.obs = adata.obs.assign(**patient_meta)
    adata.X = scipy.sparse.csr_matrix(adata.X)
    return adata


# %%
next(full_meta.iterrows())

# %%
adatas = process_map(load_counts, [r for _, r in full_meta.iterrows()], max_workers=16)

# %%
for adata in adatas:
    print(adata.shape)

# %%
for adata in adatas: 
    adata.obs_names = adata.obs["sample_id"] + "_" + adata.obs_names

# %%
adata = anndata.concat(adatas, join="outer")

# %%
adata

# %%
assert adata.obs_names.is_unique

# %%
adata.obs

# %%
adata.var

# %%
adata.write_h5ad("../data/02_anndata/nmp_liver_raw_counts_unfiltered.h5ad", compression="lzf")

# %%
adata.X

# %%
adata.X.data[:10]

# %%
