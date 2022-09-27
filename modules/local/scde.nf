process DE_DESEQ2 {
    /**
     * For standard pseudobulk analysis
     */
    cpus 2
    container "${projectDir}/containers/deseq2.sif"

    input:
    tuple val(id), path(counts), path(samplesheet)
    val(comparison)
    val(condition_col)
    val(covariate_formula)

    output:
    path("*_DESeq2_result.tsv"), emit: de_res

    script:
    def comparison_params = (comparison == "sum2zero") ? "--sum2zero" : """--c1 "${comparison[0]}" --c2 "${comparison[1]}" """
    """
    # nextflow takes care of escpaping filenames already. Do not quote.
    run_deseq2.R $counts $samplesheet \\
        --cond_col "$condition_col" \\
        $comparison_params \\
        --resDir "." \\
        --prefix "$id" \\
        --covariate_formula "$covariate_formula" \\
        --cpus $task.cpus > "${id}_deseq2.log"
    """
}
