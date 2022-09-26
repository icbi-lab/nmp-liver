"""
Functions to find and test gene signatures.

Inspired by MCPCounter and BioQC.

The original cutoffs used in the MCP counter paper are
 * fold change >= ``min_fc`` (default=2)
 * specific fold change >= ``min_sfc`` (default=1.5)
 * AUC ROC >= ``min_auc`` (default=0.97)
"""

from typing import Callable, Dict, List, Mapping, Optional, Sequence
import numpy as np
from sklearn.metrics import roc_auc_score
import sklearn.model_selection
from anndata import AnnData
import itertools
import scipy.stats
from .pseudobulk import pseudobulk
import scanpy as sc
import altair as alt
import pandas as pd
from tqdm.contrib.concurrent import process_map


def fold_change(
    adata, *, obs_col, positive_class, key_added="fold_change", inplace=True
):
    """
    Compute the fold change between the
    positive and negative samples for a single gene `G`.

    According to `Becht et al.`_ the fold change is defined as
    .. math::
        FC = X - \overline{X}
    where :math:`X` is the mean of positive and :math:`\overline{X}` is the mean of the negative samples.

    Parameters
    ----------
    adata
        annotated data matrix. Values in X are expected to be log-transformed and normalized.
    obs_col
        Column in adata.obs with the annotation to compare
    positive_class
        Value in `adata.obs[obs_col]` to be considered the positive class
    key_added
        Key under which the result will be stored in `adata.var`
    inplace
        If True, store the result in adata.var. Otherwise return result.

    Returns
    -------
    fold change

    .. _Becht et al.:
        http://dx.doi.org/10.1186/s13059-016-1070-5
    """
    positive_mask = adata.obs[obs_col] == positive_class
    fold_change = np.mean(adata.X[positive_mask, :], axis=0) - np.mean(
        adata.X[~positive_mask, :], axis=0
    )

    if inplace:
        adata.var[key_added] = fold_change
    else:
        return fold_change


def specific_fold_change(
    adata, *, obs_col, positive_class, key_added="specific_fold_change", inplace=True
):
    """
    Compute the specific fold change of the positive class with respect to all other classes
    for a single gene `G`.
    According to `Becht et al.`_ the specific fold change is defined as
    .. math::
        sFC = (X - \overline{X}_{min})/(\overline{X}_{max} - \overline{X}_{min})
    where :math:`X` is the mean of positive and :math:`\overline{X}_{max}` is the maximum mean over
    all negative classes and :math:`\overline{X}_{min}` is the minimal mean over all negative classes.

    Parameters
    ----------
    adata
        annotated data matrix. Values in X are expected to be log-transformed and normalized
    obs_col
        Column in adata.obs with the annotation to compare
    positive_class
        Value in `adata.obs[obs_col]` to be considered the positive class
    key_added
        Key under which the result will be stored in `adata.var`
    inplace
        If True, store the result in adata.var. Otherwise return result.

    Returns
    -------
    specific fold change

    .. _Becht et al.:
        http://dx.doi.org/10.1186/s13059-016-1070-5
    """
    neg_classes = [x for x in adata.obs[obs_col].unique() if x != positive_class]
    mean_per_class = np.hstack(
        [
            np.mean(adata.X[adata.obs[obs_col] == tmp_class, :], axis=0)[:, np.newaxis]
            for tmp_class in neg_classes
        ]
    )
    x_min = np.min(mean_per_class, axis=1)
    x_max = np.max(mean_per_class, axis=1)
    spec_fc = np.divide(
        np.mean(adata.X[adata.obs[obs_col] == positive_class, :], axis=0) - x_min,
        x_max - x_min,
    )

    if inplace:
        adata.var[key_added] = spec_fc
    else:
        return spec_fc


def roc_auc(
    adata,
    *,
    obs_col,
    positive_class,
    key_added="roc_auc",
    inplace=True,
):
    """
    Compute the Receiver Operator Characteristics Area under the Curve (ROC AUC) for a single gene `G`.
    This tells how well the gene discriminates between the two classes.
    This is a wrapper for the scikit-learn `roc_auc_score`_ function.

    Parameters
    ----------
    adata
        annotated data matrix. Values in X are expected to be log-transformed and normalized
    obs_col
        Column in adata.obs with the annotation to compare
    positive_class
        Value in `adata.obs[obs_col]` to be considered the positive class
    key_added
        Key under which the result will be stored in `adata.var`
    inplace
        If True, store the result in adata.var. Otherwise return result.

    Returns
    -------
    roc auc score

    .. _roc_auc_score:
        http://scikit-learn.org/stable/modules/generated/sklearn.metrics.roc_auc_score.html#sklearn.metrics.roc_auc_score
    """

    def _roc_auc(arr, mask):
        return roc_auc_score(mask, arr)

    positive_mask = adata.obs[obs_col] == positive_class
    roc_auc = np.apply_along_axis(_roc_auc, 0, adata.X, positive_mask)

    if inplace:
        adata.var[key_added] = roc_auc
    else:
        return roc_auc


def plot_metric_strip(
    pb, markers, *, metric_key="auroc", metric_title="AUROC", domain=[0.5, 1], top=10
):
    """A one-row heatmap accompaniying the heatmap generated by plot_markers showing the marker
    quality."""
    metric = []
    genes = []
    cts = []
    for ct, tmp_markers in markers.items():
        for gene in tmp_markers[:top]:
            metric.append(pb.var.loc[gene, f"{ct}_{metric_key}"])
            genes.append(gene)
            cts.append(ct)

    tmp_df = pd.DataFrame().assign(marker=genes, cell_type=cts, metric=metric)

    return (
        alt.Chart(tmp_df)
        .mark_rect()
        .encode(
            x=alt.X("marker", sort=None),
            color=alt.Color(
                "metric",
                scale=alt.Scale(scheme="viridis", domain=domain),
                title=metric_title,
            ),
        )
    )


def plot_markers(pb, groupby, markers, *, top=10, cmap="bwr", **kwargs):
    return sc.pl.matrixplot(
        pb,
        groupby=groupby,
        var_names={k: v[:top] for k, v in markers.items()},
        cmap=cmap,
        **kwargs,
    )
