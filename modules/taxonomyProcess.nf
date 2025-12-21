#!/usr/bin/env nextflow

// Process to assign taxonomy using QIIME2
process taxonomyProcess {

    conda 'envs/qiime2-env.yml'

    publishDir "${params.out_dir}/taxonomy", mode: 'copy'

    input:
    tuple val(srr_id), file(trimmed_r1), file(trimmed_r2)

    output:
    tuple val(srr_id), 
    file("${srr_id}_feature_table.tsv"), 
    file("${srr_id}_taxonomy.tsv")

    script:
    """
    touch ${srr_id}_manifest.tsv

    echo "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > ${srr_id}_manifest.tsv
    echo "${srr_id}\t${params.pwd}/${params.out_dir}/trimmed/${trimmed_r1}\t${params.pwd}/${params.out_dir}/trimmed/${trimmed_r2}" >> ${srr_id}_manifest.tsv

    qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path ${srr_id}_manifest.tsv \
        --input-format PairedEndFastqManifestPhred33V2 \
        --output-path ${srr_id}_demux.qza

    qiime dada2 denoise-paired \
        --i-demultiplexed-seqs ${srr_id}_demux.qza \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --o-table ${srr_id}_feature_table.qza \
        --o-representative-sequences ${srr_id}_rep_seqs.qza \
        --o-denoising-stats ${srr_id}_denoising_stats.qza

    qiime feature-classifier classify-sklearn \
        --i-classifier ${params.qiime_classifier} \
        --i-reads ${srr_id}_rep_seqs.qza \
        --o-classification ${srr_id}_taxonomy.qza

    qiime tools export \
        --input-path ${srr_id}_taxonomy.qza \
        --output-path ${srr_id}_taxonomy_exported

    cp ${srr_id}_taxonomy_exported/taxonomy.tsv ${srr_id}_taxonomy.tsv

    qiime tools export \
        --input-path ${srr_id}_feature_table.qza \
        --output-path ${srr_id}_feature_table_exported

    biom convert -i ${srr_id}_feature_table_exported/feature-table.biom -o ${srr_id}_feature_table.tsv --to-tsv
    """
}