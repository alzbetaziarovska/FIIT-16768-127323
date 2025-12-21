#!/usr/bin/env nextflow

// Process to create FastQC report for each file
process fastqcProcess {

    conda 'envs/core-env.yml'

    publishDir "${params.out_dir}/fastqc_${state}", mode: 'copy'

    input:
    tuple val(srr_id), file(fastq_file1), file(fastq_file2)
    val state

    output:
    tuple file("${srr_id}_1_${state}_fastqc.zip"), file("${srr_id}_1_${state}_fastqc.html"), 
          file("${srr_id}_2_${state}_fastqc.zip"), file("${srr_id}_2_${state}_fastqc.html")

    script:
    """
    mkdir -p ${params.out_dir}/fastqc_${state}
    if [ "${state}" = "raw" ]; then
        ln -s ${fastq_file1} ${srr_id}_1_raw.fastq.gz
        ln -s ${fastq_file2} ${srr_id}_2_raw.fastq.gz
        fastqc ${srr_id}_1_raw.fastq.gz
        fastqc ${srr_id}_2_raw.fastq.gz
    else
        fastqc ${fastq_file1}
        fastqc ${fastq_file2}
    fi
    """
}
