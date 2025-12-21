#!/usr/bin/env nextflow

// Process to read accession list and download data using SRA tools
process accListDownloadProcess {

    conda 'envs/core-env.yml'

    publishDir "${params.out_dir}/data", mode: 'copy'

    input:
    val srr_id

    output:
    tuple val(srr_id), file("${srr_id}_1.fastq.gz"), file("${srr_id}_2.fastq.gz")

    script:
    """
    mkdir -p ${params.out_dir}/data
    fasterq-dump --split-files ${srr_id}
    gzip ${srr_id}_1.fastq ${srr_id}_2.fastq
    """
}