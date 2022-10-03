from scanpy_helpers.cell2cell import CpdbAnalysis
import pytest
from anndata import AnnData
import pandas as pd
from pathlib import Path
from . import TESTDATA_PATH


@pytest.fixture
def adata_cell_types():
    df = pd.DataFrame.from_records(
        [
            ["T CD8", "P1", 10, 0, 0, 0, 0],
            ["T CD8", "P1", 15, 0, 0, 0, 0],
            ["T CD8", "P1", 20, 0, 0, 0, 0],
            ["T CD8", "P2", 6, 2, 0, 0, 0],
            ["T CD8", "P2", 4, 4, 0, 0, 0],
            ["receptor_cells", "P1", 0, 20, 30, 40, 0],
            ["receptor_cells", "P1", 0, 25, 35, 45, 0],
            ["receptor_cells", "P1", 0, 15, 25, 35, 0],
        ],
        columns=["cell_type_col", "patient", "CD8A", "PDGFRB", "ITGA8", "ITGB1", "FN1"],
    )
    return AnnData(
        df.drop(columns=["cell_type_col", "patient"]),
        obs=df.loc[:, ["cell_type_col", "patient"]],
    )


@pytest.fixture
def cellchatdb():
    return pd.read_csv(TESTDATA_PATH / "cellchatdb_minimal.tsv", sep="\t", comment="#")


@pytest.fixture
def cpdba(adata_cell_types, cellchatdb):
    return CpdbAnalysis(
        cellchatdb,
        adata_cell_types,
        pseudobulk_group_by=["patient"],
        cell_type_column="cell_type_col",
        min_obs=1,
    )


def test_find_expressed_genes(cpdba):
    """Check that computing the mean expression and mean fraction of cells expressing a gene works as expected"""
    assert cpdba.expressed_genes.to_records(index=False).tolist() == [('T CD8', 'CD8A', 1.0, pytest.approx(13.5805)), ('receptor_cells', 'CD8A', 0.0, 0.0), ('T CD8', 'PDGFRB', 0.5, 6.417342185974121), ('receptor_cells', 'PDGFRB', 1.0, 12.311437606811523), ('T CD8', 'ITGA8', 0.0, 0.0), ('receptor_cells', 'ITGA8', 1.0, 12.716901779174805), ('T CD8', 'ITGB1', 0.0, 0.0), ('receptor_cells', 'ITGB1', 1.0, 13.004582405090332), ('T CD8', 'FN1', 0.0, 0.0), ('receptor_cells', 'FN1', 0.0, 0.0)]


@pytest.mark.parametrize("de_genes,expected_interactions", [
    [], []
])
def test_find_significant_interactions(cpdba, de_genes):
