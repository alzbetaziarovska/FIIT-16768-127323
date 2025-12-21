#!/usr/bin/env nextflow

// Process to trim adapters and low-quality bases from FASTQ files using Trimmomatic
process trimProcess {

    conda 'envs/core-env.yml'

    publishDir "${params.out_dir}/trimmed", mode: 'copy'

    input:
    tuple val(srr_id), file(fastq_file1), file(fastq_file2)

    output:
    tuple val(srr_id), file("${srr_id}_1_trimmed.fastq.gz"), file("${srr_id}_2_trimmed.fastq.gz")

    script:
    """
    mkdir -p ${params.out_dir}/trimmed
    cutadapt -a ^${params.adapter} -o ${srr_id}_1_trimmed.fastq ${fastq_file1}
    cutadapt -a ^${params.adapter} -o ${srr_id}_2_trimmed.fastq ${fastq_file2}
    gzip ${srr_id}_1_trimmed.fastq ${srr_id}_2_trimmed.fastq
    """
}