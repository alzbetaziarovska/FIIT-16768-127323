#!/usr/bin/env nextflow

// Process to perform statistical analysis
process statisticsProcess {
    maxForks 1
    conda 'envs/statistics-env.yml'
    publishDir "${params.out_dir}/statistics", mode: 'copy'

    input:
    val tuple_list

    output:
    file "merged_abundance.csv"
    file "merged_taxonomy.csv"
    file "${params.name}_*.pdf"

    script:
    // Collect table files and sample IDs as comma-separated strings
    def table_files_input = tuple_list.collect { it[1] }.join(",")
    def sample_ids_input = tuple_list.collect { it[0] }.join(",")

    """
    cat <<'EOF' > temp_script.R
    library(dplyr)
    library(vegan)
    library(ggplot2)
    library(phyloseq)
    library(stats)
    library(tidyr)
    library(ggrepel)

    table_files <- strsplit(Sys.getenv("TABLE_FILES"), ",")[[1]]
    sample_ids <- strsplit(Sys.getenv("SAMPLE_IDS"), ",")[[1]]
    top_taxa <- as.numeric(Sys.getenv("TOP_TAXA"))
    output_prefix <- Sys.getenv("OUTPUT_PREFIX")
    
    list_dfs <- list()
    taxonomy_list <- list()

    for (f in table_files) {
    df <- read.delim(f, check.names=FALSE, stringsAsFactors=FALSE)
    
    otu_col <- colnames(df)[1]
    sample_col <- colnames(df)[2]
    
    df_sub <- df %>%
        select(all_of(c(otu_col, sample_col))) %>%
        rename(OTU_ID = !!otu_col, !!sample_col := sample_col)
    
    df_sub[[sample_col]] <- as.numeric(df_sub[[sample_col]])
    list_dfs[[f]] <- df_sub
    
    tax_sub <- df %>%
        select(all_of(c(otu_col, "Taxon"))) %>%
        rename(OTU_ID = !!otu_col)
    taxonomy_list[[f]] <- tax_sub
    }

    merged_df <- Reduce(function(x, y) full_join(x, y, by="OTU_ID"), list_dfs)
    merged_df[is.na(merged_df)] <- 0

    taxonomy_df <- bind_rows(taxonomy_list) %>%
    distinct(OTU_ID, .keep_all=TRUE) %>%
    separate(Taxon,
            into=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
            sep=";", fill="right", extra="drop") %>%
    mutate(across(Kingdom:Species, ~gsub("^[kpcofgsd]__+", "", .)))

    taxonomy_df <- taxonomy_df[!grepl("Chordata", taxonomy_df[["Phylum"]], ignore.case = TRUE), ]
    taxonomy_df <- taxonomy_df[, -ncol(taxonomy_df)]

    merged_df <- merged_df[merged_df[["OTU_ID"]] %in% taxonomy_df[["OTU_ID"]], ]
    colnames(merged_df) <- c("OTU_ID", sample_ids)

    otu_ids <- merged_df[["OTU_ID"]]
    short_codes <- paste0("O", sprintf("%03d", seq_along(otu_ids)))
    otu_map <- data.frame(OTU_ID = otu_ids, ShortID = short_codes, stringsAsFactors = FALSE)

    merged_df <- merged_df %>%
    left_join(otu_map, by = "OTU_ID") %>%
    select(ShortID, everything(), -OTU_ID) %>%
    rename(OTU_ID = ShortID)

    taxonomy_df <- taxonomy_df %>%
    left_join(otu_map, by = "OTU_ID") %>%
    select(ShortID, everything(), -OTU_ID) %>%
    rename(OTU_ID = ShortID)

    write.table(merged_df, file="merged_abundance.csv", sep=";", quote=FALSE, row.names=FALSE)
    write.table(taxonomy_df, file="merged_taxonomy.csv", sep=";", quote=FALSE, row.names=FALSE)

    merged_otu_table <- read.table("merged_abundance.csv", header=TRUE, sep=";")

    rowSums(merged_df[,-1])

    otu_mat <- as.matrix(merged_otu_table[, -1])
    rownames(otu_mat) <- merged_otu_table[["OTU_ID"]]
    otu_table_obj <- otu_table(otu_mat, taxa_are_rows = TRUE)

    # Plotting
    p1 <- plot_richness(otu_table_obj, measures=c("Shannon")) +
          labs(x="Sample ID") + geom_point(color="red", size=2)

    top_taxa_names <- names(sort(taxa_sums(otu_table_obj), decreasing = TRUE))[1:top_taxa]
    otu_top <- prune_taxa(top_taxa_names, otu_table_obj)

    p2 <- plot_bar(otu_table_obj) +
          labs(x="Sample ID") +
          geom_text(aes(label = OTU), size = 3, position = position_stack(vjust = 0.5)) +
          theme(legend.position = "none")

    p3 <- plot_heatmap(otu_top, method = "NMDS", distance = "bray")

    p4 <- plot_richness(otu_table_obj, measures=c("Chao1")) +
          labs(x="Sample ID") + geom_point(color="red", size=2)

    otu_mat_t <- t(otu_mat)
    d <- vegdist(otu_mat_t, method = "bray", na.rm = TRUE)
    cluster <- hclust(d, method = "average")
    dend <- as.dendrogram(cluster)

    pca_otu <- as.data.frame(otu_table_obj)
    pca_otu_t <- t(pca_otu)
    pca_res <- rda(pca_otu_t)
    pca_scores <- scores(pca_res, display = "sites")
    pca_df <- as.data.frame(pca_scores)
    pca_df[["SampleID"]] <- rownames(pca_df)

    pdf(paste0(output_prefix, "_shannon.pdf"), width = nsamples(otu_table_obj), height = 8)
    print(p1); dev.off()

    pdf(paste0(output_prefix, "_chao1.pdf"), width = nsamples(otu_table_obj), height = 8)
    print(p4); dev.off()

    pdf(paste0(output_prefix, "_bar.pdf"), width = nsamples(otu_table_obj), height = 10)
    print(p2); dev.off()

    pdf(paste0(output_prefix, "_heatmap.pdf"), width = (top_taxa%/%2), height = (top_taxa%/%2))
    print(p3); dev.off()

    pdf(paste0(output_prefix, "_dendrogram.pdf"), width = 10, height = top_taxa)
    par(mar = c(5, 4, 4, 8) + 0.1)
    plot(dend, horiz=TRUE, main="Dendrogram of Samples", xlab="Bray-Curtis Distance", ylab="Sample ID", yaxs = "i")
    dev.off()

    pdf(paste0(output_prefix, "_pca.pdf"), width = 10, height = 8)
    ggplot(pca_df, aes(x = PC1, y = PC2)) +
        geom_point(size = 3, color = "red") +
        geom_text_repel(aes(label = SampleID), size = 3) +
        theme_minimal() +
        labs(title = "PCA plot of samples", x = "PC1", y = "PC2")
    dev.off()
    EOF

    # Export environment variables for R
    export TABLE_FILES="${table_files_input}"
    export SAMPLE_IDS="${sample_ids_input}"
    export TOP_TAXA=${params.top_taxa}
    export OUTPUT_PREFIX="${params.name}"

    # Run the R script
    Rscript temp_script.R
    """
}