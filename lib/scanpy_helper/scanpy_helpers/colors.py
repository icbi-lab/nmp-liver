import altair as alt
import pandas as pd


def set_scale_anndata(adata, column, palette=None):
    if palette is None:
        palette = column

    adata._sanitize()

    tmp_cols = getattr(COLORS, palette)
    adata.uns[f"{column}_colors"] = [
        tmp_cols[cat] for cat in adata.obs[column].cat.categories
    ]


def altair_scale(variable, *, data=None, data_col=None, **kwargs):
    """
    Discrete color scale for altair based on our color definitions.

    Parameters:
    -----------
    variable
        name of the color scale
    data
        Data frame used for the chart. If specified, will only show values that actually occur in the data.
    data_col
        If specified, check this column in `data` instead of `variable`

    Returns
    -------
    Altair color scale
    """
    tmp_colors = getattr(COLORS, variable)
    if data is not None:
        data_col = variable if data_col is None else data_col
        tmp_colors = {k: tmp_colors[k] for k in sorted(data[data_col].unique())}

    return alt.Scale(
        domain=list(tmp_colors.keys()),
        range=list(tmp_colors.values()),
        **kwargs,
    )


def altair_scale_mpl(scheme, **kwargs):
    """
    Use a continuous color scheme from mpl with altair
    """
    from matplotlib import cm
    from matplotlib.colors import to_hex

    return alt.Scale(
        range=[to_hex(x) for x in cm.get_cmap(scheme, 1000)(range(1000))], **kwargs
    )


def plot_palette(variable):
    """Display a palette"""
    tmp_cols = getattr(COLORS, variable)
    return (
        alt.Chart(
            pd.DataFrame.from_dict(tmp_cols, orient="index", columns=["color"])
            .reset_index()
            .rename(columns={"index": variable})
        )
        .mark_rect(height=40, width=30)
        .encode(
            x=alt.X(variable),
            color=alt.Color(variable, scale=altair_scale(variable), legend=None),
        )
    )


def plot_all_palettes():
    return alt.vconcat(
        *[plot_palette(v) for v in COLORS.__dict__.keys() if not v.startswith("_")]
    ).resolve_scale(color="independent")


class COLORS:
    timepoint = {
        "T0": "#edf8b1",
        "T1": "#7fcdbb",
        "T2": "#2c7fb8",
    }
    patient_id = {
        "P1": "#000000",
        "P2": "#E69F00",
        "P3": "#56B4E9",
        "P4": "#009E73",
        "P5": "#F0E442",
        "P6": "#0072B2",
        "P7": "#D55E00",
        "P8": "#CC79A7",
    }
    cell_type = {
        "B cells": "#1f77b4",
        "Cholangiocytes": "#ff7f0e",
        "Endothelial cells": "#279e68",
        "Hepatocytes": "#aa40fc",
        "Mast cells": "#8c564b",
        "NK cell": "#e377c2",
        "Neutrophils": "#b5bd61",
        "Plasma cells": "#17becf",
        "Progenitor": "#aec7e8",
        "T cell CD4": "#ffbb78",
        "T cell CD8": "#98df8a",
        "T cell CD8 NK-like": "#ff9896",
        "cDC": "#c5b0d5",
        "monocytic lineage": "#c49c94",
        "pDC": "#f7b6d2",
    }
    cell_type_coarse = {
        "B cells": "#1f77b4",
        "Cholangiocytes": "#ff7f0e",
        "Endothelial cells": "#279e68",
        "Hepatocytes": "#d62728",
        "Mast cells": "#aa40fc",
        "NK cell": "#8c564b",
        "Neutrophils": "#e377c2",
        "Plasma cells": "#b5bd61",
        "Progenitor": "#17becf",
        "T cell": "#aec7e8",
        "cDC": "#ffbb78",
        "monocytic lineage": "#98df8a",
        "pDC": "#ff9896",
    }
