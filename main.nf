#!/usr/bin/env nextflow

// Pipeline parameters
params.pwd = "$PWD" // Current working directory
params.data_file = "SRR_Acc_list.txt" // Path to the file containing SRR IDs
params.data_dir = "data" // Path to the directory where already downloaded data is stored
params.download = true // Set to true to download data from AccList, false to use existing files
params.out_dir = "results" // Output directory for results
params.name = "run1" // Name of the project, used in output file names
params.adapter = "TACGGRAGGCAGCAG...AGGGTATCTAATCCT" // Adapter sequence for trimming
params.qiime_classifier = "${params.pwd}/silva-138-99-nb-classifier.qza" // Path to pretrained QIIME2 classifier
params.top_taxa = 20 // Number of top taxa to display in the statistics report (heatmap) and double for dendrogram

// Include processes from modules
include { accListDownloadProcess } from './modules/accListDownloadProcess.nf'
include { fastqcProcess as fastqcRaw } from './modules/fastqcProcess.nf'
include { fastqcProcess as fastqcTrimmed } from './modules/fastqcProcess.nf'
include { multiqcProcess as multiqcRaw } from './modules/multiqcProcess.nf'
include { multiqcProcess as multiqcTrimmed } from './modules/multiqcProcess.nf'
include { trimProcess } from './modules/trimProcess.nf'
include { taxonomyProcess } from './modules/taxonomyProcess.nf'
include { statisticsProcess } from './modules/statisticsProcess.nf'
include { convertProcess } from './modules/convertProcess.nf'


// Define the workflow
workflow {

    if (params.download) {
        srr_ids = Channel
            .fromPath(params.data_file)
            .splitText()
            .map { it.trim() }

        accListDownloadProcess(srr_ids)
    } else {
        // load existing .fastq.gz files from params.data_dir
        // Expecting files named as <SRR_ID>_1.fastq.gz and <SRR_ID>_2.fastq.gz for paired-end reads

        srr_files = Channel
            .fromPath("${params.data_dir}/*_1.fastq.gz")
            .map { file ->
                def srr_id = file.getSimpleName().replaceAll(/_1$/, "")
                tuple(srr_id, file, file.getParent().resolve("${srr_id}_2.fastq.gz"))
            }
        
    }

    input_files = params.download ? accListDownloadProcess.out : srr_files

    // Run FastQC on raw files
    fastqcRaw(input_files, "raw")

    // Run Trimmomatic on the downloaded files
    trimmed_reads = trimProcess(input_files)
    fastqcTrimmed(trimmed_reads, "trimmed")

    // Collect FastQC output files and run MultiQC
    fastqcFilesRaw = fastqcRaw.out
        .collect()
        .flatten()
        .collect()

    fastqcFilesTrimmed = fastqcTrimmed.out
        .collect()
        .flatten()
        .collect()

    // Compile MultiQC report
    multiqcRaw(fastqcFilesRaw, "raw")
    multiqcTrimmed(fastqcFilesTrimmed, "trimmed")

    // Run taxonomy analysis on the processed files
    taxonomyProcess(trimmed_reads)

    // Convert QIIME2 artifacts to tsv files and collect them
    convertProcess(taxonomyProcess.out)
    tables = convertProcess.out
       .toList()

    // Generate statistics report
    statisticsProcess(tables)

}

