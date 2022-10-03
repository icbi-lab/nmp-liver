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
            ["T CD8", "P1", 10, 0, 0, 0, 0, 5],
            ["T CD8", "P1", 15, 0, 1, 0, 0, 0],
            ["T CD8", "P1", 20, 0, 0, 0, 0, 0],
            ["T CD8", "P2", 6, 2, 1, 0, 0, 0],
            ["T CD8", "P2", 4, 4, 0, 0, 0, 0],
            ["receptor_cells", "P1", 0, 20, 30, 40, 1, 0],
            ["receptor_cells", "P1", 0, 25, 35, 45, 0, 0],
            ["receptor_cells", "P1", 0, 15, 25, 35, 0, 0],
        ],
        columns=[
            "cell_type_col",
            "patient",
            "CD8A",
            "PDGFRB",
            "ITGA8",
            "ITGB1",
            "FN1",
            "SDC4",
        ],
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
    # The mean expression is a geometric mean (i.e. mean of log-transformed values)
    assert cpdba.expressed_genes.to_records(index=False).tolist() == [
        # cell_type, gene, fraction_expressed, mean_expressed
        ("T CD8", "CD8A", 1.0, pytest.approx(13.487, abs=0.001)),
        ("receptor_cells", "CD8A", 0.0, pytest.approx(0.0, abs=0.001)),
        ("T CD8", "PDGFRB", 0.5, pytest.approx(6.3870, abs=0.001)),
        ("receptor_cells", "PDGFRB", 1.0, pytest.approx(12.307, abs=0.001)),
        (
            "T CD8",
            "ITGA8",
            pytest.approx(0.4166, abs=0.001),
            pytest.approx(10.433, abs=0.001),
        ),
        ("receptor_cells", "ITGA8", 1.0, pytest.approx(12.713, abs=0.001)),
        ("T CD8", "ITGB1", 0.0, pytest.approx(0.0, abs=0.001)),
        ("receptor_cells", "ITGB1", 1.0, pytest.approx(13.000, abs=0.001)),
        ("T CD8", "FN1", 0.0, pytest.approx(0.0, abs=0.001)),
        (
            "receptor_cells",
            "FN1",
            pytest.approx(0.333, abs=0.001),
            pytest.approx(8.214, abs=0.001),
        ),
        (
            "T CD8",
            "SDC4",
            pytest.approx(0.1666, abs=0.001),
            pytest.approx(5.74656, abs=0.001),
        ),
        ("receptor_cells", "SDC4", 0.0, 0.0),
    ]


@pytest.mark.parametrize(
    "de_genes,params,expected_interactions",
    [
        (
            ["PDGFB", "doesntexist"],
            {},
            [("T CD8", "PDGFB", "PDGFRB"), ("receptor_cells", "PDGFB", "PDGFRB")],
        ),
        (
            ["FN1"],
            {},
            [("T CD8", "FN1", "SDC4")],
        ),
        (
            ["FN1"],
            {"max_pvalue": 0.000001},
            [],
        ),
        (
            ["FN1"],
            {"min_abs_fc": 1.1},
            [],
        ),
        (
            ["FN1"],
            {"min_frac_expressed": 0.5},
            [],
        ),
        (
            ["VTN"],
            {"complex_policy": "ignore"},
            [],
        ),
        (
            ["SDC4", "doesntexist", "CD8A", "PDGFB"],
            {"de_genes_mode": "receptor"},
            [("receptor_cells", "FN1", "SDC4")],
        ),
        (
            ["VTN"],
            {"complex_policy": "explode"},
            [
                ("T CD8", "VTN", "ITGA8"),
                ("receptor_cells", "VTN", "ITGA8"),
                ("receptor_cells", "VTN", "ITGB1"),
            ],
        ),
    ],
)
def test_find_significant_interactions(cpdba, de_genes, params, expected_interactions):
    # Don't care about pvalue here for now
    de_res = pd.DataFrame().assign(gene_symbol=de_genes, log_fc=1, p_value=0.0001)
    interactions = cpdba.significant_interactions(
        de_res,
        pvalue_col="p_value",
        fc_col="log_fc",
        gene_symbol_col="gene_symbol",
        **params
    )
    assert (
        sorted(
            interactions.loc[
                :, ["cell_type_col", "source_genesymbol", "target_genesymbol"]
            ]
            .to_records(index=False)
            .tolist()
        )
        == expected_interactions
    )
