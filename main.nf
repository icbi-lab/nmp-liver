#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { check_samplesheet } from './modules/local/check_samplesheet'
include { SCQC } from "./modules/local/scqc/main"

workflow {

    ch_samples = Channel.from(check_samplesheet("${projectDir}/tables/samplesheet_scrnaseq_qc.csv", projectDir))

    SCQC(
        [
            file("${projectDir}/modules/local/scqc/scqc-notebook.py", checkIfExists: true),
            file("${projectDir}/modules/local/scqc/qc_plots.py", checkIfExists: true)
        ],
        ch_samples
    )

}
