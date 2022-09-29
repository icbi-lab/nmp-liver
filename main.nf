#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { check_samplesheet } from './modules/local/check_samplesheet'
include { SCQC } from "./modules/local/scqc/main"
include { JUPYTERNOTEBOOK as JUPYTER_SCVI } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_CELL_TYPES } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_MYELOID } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_NEUTRO } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_DE_ANALYSIS } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_OVERVIEW_PLOTS } from "./modules/local/jupyternotebook/main"
include { DE_DESEQ2 as DESEQ_T0_T1 } from "./modules/local/scde.nf"
include { DE_DESEQ2 as DESEQ_LT } from "./modules/local/scde.nf"
include { DE_DESEQ2 as DESEQ_ECD } from "./modules/local/scde.nf"

workflow {

    ch_samples = Channel.from(check_samplesheet("${projectDir}/tables/samplesheet_scrnaseq_qc.csv", projectDir))

    SCQC(
        [
            file("${projectDir}/modules/local/scqc/scqc-notebook.py", checkIfExists: true),
            file("${projectDir}/modules/local/scqc/qc_plots.py", checkIfExists: true)
        ],
        ch_samples
    )

    ch_adata_qc = SCQC.out.adata.map {id, adata -> adata }
    JUPYTER_SCVI(
        Channel.value([
            [id: "03_scvi"],
            file("${projectDir}/analyses/03_scvi.py", checkIfExists: true)
        ]),
        ch_adata_qc.map { adata -> ["adata_path": adata.name] },
        ch_adata_qc
    )

    ch_cell_type_annotation = JUPYTER_SCVI.out.artifacts.flatten().filter{
        it -> it.name.contains("adata_scvi_doublet_filtered")
    }.concat(
        Channel.fromPath("${projectDir}/tables/cell_type_markers.csv")
    ).collect()
    JUPYTER_CELL_TYPES(
        Channel.value([
            [id: "04_cell_type_annotation"],
            file("${projectDir}/analyses/04_cell_type_annotation.py", checkIfExists: true)
        ]),
        ch_cell_type_annotation.map{ adata, markers -> [
            "adata_path": adata.name,
            "marker_genes_path": markers.name
        ]},
        ch_cell_type_annotation
    )

    ch_adata_cell_types = JUPYTER_CELL_TYPES.out.artifacts.flatten().filter{
        it -> it.name.contains("adata_cell_types")
    }

    JUPYTER_DE_ANALYSIS(
        Channel.value([
            [id: "20_de_analysis"],
            file("${projectDir}/analyses/20_de_analysis.py", checkIfExists: true)
        ]),
        ch_adata_cell_types.map{ adata -> [
            "adata_path": adata.name
        ]},
        ch_adata_cell_types
    )

    ch_deseq = JUPYTER_DE_ANALYSIS.out.artifacts.flatten().map{
        it -> [it.name.split("\\.")[0], it]
    }.groupTuple(sort: true).map{
        id, files -> [id, files[0], files[1]]
    }

    DESEQ_T0_T1(
        ch_deseq.filter{ it -> it[0].contains("t0_vs_t1") },
        ["T1", "T0"],
        "timepoint",
        "+ patient_id"
    )
    DESEQ_LT(
        ch_deseq.filter{ it -> it[0].contains("LT_")},
        ["Yes", "No"],
        "LT",
        ""
    )
    DESEQ_ECD(
        ch_deseq.filter{ it -> it[0].contains("ECD_")},
        ["Yes", "No"],
        "ECD",
        ""
    )

    ch_overview_plots = ch_adata_cell_types.concat(
        Channel.fromPath("${projectDir}/tables/cell_type_markers.csv")
    ).collect()
    JUPYTER_OVERVIEW_PLOTS(
        Channel.value([
            [id: "90_overview_plots"],
            file("${projectDir}/analyses/90_overview_plots.py", checkIfExists: true)
        ]),
        ch_overview_plots.map{ adata, markers -> [
            "adata_path": adata.name,
            "marker_genes_path": markers.name
        ]},
        ch_overview_plots
    )

}
