#!/usr/bin/env nextflow

// Process to convert feature and taxonomy tables into a unified format
process convertProcess {
    
    conda 'envs/statistics-env.yml'

    publishDir "${params.out_dir}/taxonomy", mode: 'copy'

    input:
    tuple val(srr_id), val(feature_table), val(taxonomy_table)

    output:
    tuple val(srr_id), file("${srr_id}_table.tsv")

    script:
    """
    cat <<'EOF' > temp_script.R

        tax_table <- Sys.getenv("TAX_TABLE")
        feature_table <- Sys.getenv("FEATURE_TABLE")
        srr_id <- Sys.getenv("SRR_ID")

        taxonomy <- read.delim(tax_table, stringsAsFactors=FALSE, check.names=FALSE)
        feature <- read.delim(feature_table, stringsAsFactors=FALSE, check.names=FALSE, skip=1)

        feature_col <- feature[[2]]
        feature_name <- colnames(feature)[2]

        taxonomy <- cbind(taxonomy[1], feature_col, taxonomy[-1])
        colnames(taxonomy)[2] <- feature_name

        colnames(taxonomy)[1] <- "OTU_ID"

        taxonomy <- taxonomy[, !(colnames(taxonomy) %in% "Confidence")]

        write.table(taxonomy, file=paste0(srr_id, "_table.tsv"), sep="\t", row.names=FALSE, quote=FALSE)
    EOF

    export TAX_TABLE="${taxonomy_table}"
    export FEATURE_TABLE="${feature_table}"
    export SRR_ID="${srr_id}"

    Rscript temp_script.R
    """

}