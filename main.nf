#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { check_samplesheet } from './modules/local/check_samplesheet'
include { SCQC } from "./modules/local/scqc/main"
include { JUPYTERNOTEBOOK as JUPYTER_SCVI } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_CELL_TYPES } from "./modules/local/jupyternotebook/main"

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

}
