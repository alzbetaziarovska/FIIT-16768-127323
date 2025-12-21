#!/usr/bin/env nextflow

// Process to generate a MultiQC report from FastQC outputs
process multiqcProcess {

    conda 'envs/multiqc-env.yml'

    publishDir "${params.out_dir}", mode: 'copy'

    input:
    path fastqc_files, stageAs: 'fastqc/*'
    val state

    output:
    file "multiqc_${state}_report.html"

    script:
    """
    multiqc . --filename multiqc_${state}_report.html
    """
}