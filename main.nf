#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { check_samplesheet } from './modules/local/check_samplesheet'
include { SCQC } from "./modules/local/scqc/main"
include { JUPYTERNOTEBOOK as JUPYTER_SCVI } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_CELL_TYPES } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_MYELOID } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_NEUTRO } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_NKT } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_DE_ANALYSIS } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_OVERVIEW_PLOTS } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_T0T1 } from "./modules/local/jupyternotebook/main"
include { JUPYTERNOTEBOOK as JUPYTER_LIVER_QUALITY } from "./modules/local/jupyternotebook/main"
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

    JUPYTER_MYELOID(
        Channel.value([
            [id: "10_myeloid_analysis"],
            file("${projectDir}/analyses/10_myeloid_analysis.py", checkIfExists: true)
        ]),
        ch_adata_cell_types.map{ adata -> [
            "adata_path": adata.name
        ]},
        ch_adata_cell_types
    )
    JUPYTER_NEUTRO(
        Channel.value([
            [id: "11_neutro_analysis"],
            file("${projectDir}/analyses/11_neutro_analysis.py", checkIfExists: true)
        ]),
        ch_adata_cell_types.map{ adata -> [
            "adata_path": adata.name
        ]},
        ch_adata_cell_types
    )
    JUPYTER_NKT(
        Channel.value([
            [id: "12_nk_and_t_analysis"],
            file("${projectDir}/analyses/12_nk_and_t_analysis.py", checkIfExists: true)
        ]),
        ch_adata_cell_types.map{ adata -> [
            "adata_path": adata.name
        ]},
        ch_adata_cell_types
    )

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

    ch_adata_m = JUPYTER_MYELOID.out.artifacts.flatten().filter{
        it -> it.name.contains("adata_m")
    }
    ch_adata_n = JUPYTER_NEUTRO.out.artifacts.flatten().filter{
        it -> it.name.contains("adata_n")
    }
    ch_adata_nkt = JUPYTER_NKT.out.artifacts.flatten().filter{
        it -> it.name.contains("adata_nkt")
    }

    ch_t0t1 = ch_adata_cell_types.concat(
        ch_adata_m,
        ch_adata_n,
        ch_adata_nkt,
        Channel.fromPath("${projectDir}/tables/dorothea_human_AB_2022-09-28.csv"),
        Channel.fromPath("${projectDir}/tables/cellchatdb_2022-09-29.tsv"),
        Channel.fromPath("${projectDir}/tables/gene_sets_hallmarks_msigdb.csv"),
        Channel.fromPath("${projectDir}/tables/gene_sets_interleukins_chemokines.xlsx")
    ).collect()
    JUPYTER_T0T1(
       Channel.value([
            [id: "40_t0_vs_t1"],
            file("${projectDir}/analyses/40_t0_vs_t1.py", checkIfExists: true)
        ]),
        ch_t0t1.map{ adata, adata_m, adata_n, adata_nkt, dorothea, cellchatdb, msigdb, gene_set_il -> [
            "adata_path": adata.name,
            "adata_myeloid_path": adata_m.name,
            "adata_neutro_path": adata_n.name,
            "adata_nkt_path": adata_nkt.name,
            "de_res_dir": ".",
            "dorothea": dorothea.name,
            "cellchatdb": cellchatdb.name,
            "msigdb": msigdb.name,
            "gene_set_il_path": gene_set_il.name
        ]},
        ch_t0t1.mix(DESEQ_T0_T1.out.de_res).collect()
    )

    ch_de_res_quality = DESEQ_ECD.out.de_res.mix(DESEQ_LT.out.de_res).collect()
    ch_liver_quality = ch_adata_cell_types.concat(
        Channel.fromPath("${projectDir}/tables/gene_sets_interleukins_chemokines.xlsx")
    ).collect()
    JUPYTER_LIVER_QUALITY(
        Channel.value([
            [id: "41_liver_quality"],
            file("${projectDir}/analyses/41_liver_quality.py", checkIfExists: true)
        ]),
        ch_liver_quality.map{ adata, gene_set_il -> [
            "adata_path": adata.name,
            "gene_set_il_path": gene_set_il.name,
            "de_res_dir": "."
        ]},
        ch_liver_quality.mix(ch_de_res_quality).collect()
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
