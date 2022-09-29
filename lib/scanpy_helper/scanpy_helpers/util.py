from tqdm.auto import tqdm
import contextlib
import os
import statsmodels.stats.multitest
import numpy as np
from anndata import AnnData
from matplotlib.patches import PathPatch
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
    for cat in categories:
        yield cat, adata[adata.obs[groupby] == cat, :].copy()


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


def adjust_box_widths(g, fac):
    """
    Adjust the withs of a seaborn-generated boxplot.

    from https://stackoverflow.com/questions/56838187/how-to-create-spacing-between-same-subgroup-in-seaborn-boxplot
    """

    # iterating through Axes instances
    for ax in g.axes:

        # iterating through axes artists:
        for c in ax.get_children():

            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5 * (xmin + xmax)
                xhalf = 0.5 * (xmax - xmin)

                # setting new width of box
                xmin_new = xmid - fac * xhalf
                xmax_new = xmid + fac * xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new

                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])
