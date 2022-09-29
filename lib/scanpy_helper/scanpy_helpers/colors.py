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
        "T0": "#1b9e77",
        "T1": "#d95f02",
        "T2": "#7570b3",
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
        "B cells": "#b5bd61",
        "Cholangiocytes": "#f7b6d2",
        "Endothelial cells": "#8c564b",
        "Hepatocytes": "#17becf",
        "Mast cells": "#ffbb78",
        "NK cell": "#aa40fc",
        "Neutrophils": "#1f77b4",
        "Plasma cells": "#e377c2",
        "Progenitor": "#333333",
        "T cell CD4": "#279e68",
        "T cell CD8": "#006d2c",
        "T cell CD8 NK-like": "#98df8a",
        "cDC": "#aec7e8",
        "monocytic lineage": "#ff7f0e",
        "pDC": "#c49c94",
    }
    cell_type_coarse = {
        "B cells": "#b5bd61",
        "Cholangiocytes": "#f7b6d2",
        "Endothelial cells": "#8c564b",
        "Hepatocytes": "#17becf",
        "Mast cells": "#ffbb78",
        "NK cell": "#aa40fc",
        "Neutrophils": "#1f77b4",
        "Plasma cells": "#e377c2",
        "Progenitor": "#333333",
        "T cell": "#279e68",
        "cDC": "#aec7e8",
        "monocytic lineage": "#ff7f0e",
        "pDC": "#c49c94",
    }
    LT = {
        "yes": "#ff7f00",
        "no": "#4daf4a",
    }
    ECD = {
        "yes": "#d95f02",
        "no": "#7570b3",
    }
