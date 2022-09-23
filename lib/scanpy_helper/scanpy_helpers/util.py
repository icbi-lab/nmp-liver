from tqdm.auto import tqdm
import contextlib
import os
import statsmodels.stats.multitest
import numpy as np
from anndata import AnnData
import scipy.sparse



def fdr_correction(df, pvalue_col="pvalue", *, key_added="fdr", inplace=False):
    """Adjust p-values in a data frame with test results using FDR correction."""
    if not inplace:
        df = df.copy()

    df[key_added] = statsmodels.stats.multitest.fdrcorrection(df[pvalue_col].values)[1]

    if not inplace:
        return df


def split_anndata(adata, groupby):
    """Split an anndata object into a dict of anndata objects based on a column in obs"""
    categories = adata.obs[groupby].unique()
    return {cat: adata[adata.obs[groupby] == cat, :].copy() for cat in tqdm(categories)}


def chunk_adatas(ad, chunksize=200):
    """Generate chunks of adata objects (by variable)"""
    for i in range(0, ad.shape[1], chunksize):
        yield ad[:, i : i + chunksize].copy()


def suppress_stdout(func):
    """Decorator to suppress stdout"""

    def wrapper(*a, **ka):
        with open(os.devnull, "w") as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*a, **ka)

    return wrapper


def _choose_mtx_rep(adata, use_raw=False, layer=None):
    is_layer = layer is not None
    if use_raw and is_layer:
        raise ValueError(
            "Cannot use expression from both layer and raw. You provided:"
            f"'use_raw={use_raw}' and 'layer={layer}'"
        )
    if is_layer:
        return adata.layers[layer]
    elif use_raw:
        return adata.raw.X
    else:
        return adata.X





def scale_range(a):
    """Scale between -1 and 1, centered around the original 0"""
    return a / max(np.abs(np.max(a)), np.abs(np.min(a)))


def scale_01(a):
    """Scale between 0 and 1"""
    return (a - np.min(a)) / np.ptp(a)
