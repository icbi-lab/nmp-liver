#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { check_samplesheet } from './modules/local/check_samplesheet'
include { SCQC } from "./modules/local/scqc/main"
include { JUPYTERNOTEBOOK as JUPYTER_SCVI } from "./modules/local/jupyternotebook/main"

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
            [id: "02_scvi"],
            file("${projectDir}/analyses/02_scvi.py", checkIfExists: true)
        ]),
        ch_adata_qc.map { adata -> ["adata_path": adata.name] },
        ch_adata_qc
    )

}
