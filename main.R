library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))

make_se <- function(counts_csv, metafile_csv, selected_times) {
  # Read in counts and metadata
  counts <- read.table(counts_csv, header = TRUE, sep = "\t", row.names = 1)
  metadata <- read.csv(metafile_csv)
  
  # Filter metadata to selected timepoints
  metadata_filtered <- metadata[metadata$timepoint %in% selected_times, ]
  
  # Set vP0 as reference factor level for timepoint
  metadata_filtered$timepoint <- factor(metadata_filtered$timepoint)
  metadata_filtered$timepoint <- relevel(metadata_filtered$timepoint, ref = "vP0")
  
  # Subset counts to only include samples in filtered metadata
  counts_filtered <- counts[, colnames(counts) %in% metadata_filtered$samplename]
  
  # Ensure column order of counts matches row order of metadata
  counts_filtered <- counts_filtered[, metadata_filtered$samplename]
  
  # Build the SummarizedExperiment object
  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts_filtered)),
    colData = DataFrame(
      samplename = metadata_filtered$samplename,
      timepoint  = metadata_filtered$timepoint
    )
  )
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)

return_deseq_res <- function(se, design) {
  # Create DESeqDataSet from SummarizedExperiment
  dds <- DESeqDataSet(se, design = design)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results as a dataframe
  res <- results(dds)
  res_df <- as.data.frame(res)
  
  # Return named list with dds object and results dataframe
  return(list(
    dds = dds,
    results = res_df
  ))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
    #' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  labeled <- deseq2_res %>%
    as.data.frame() %>%                        # <-- ensures rownames are preserved
    tibble::rownames_to_column(var = "genes") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(volc_plot_status = dplyr::case_when(
      padj < padj_threshold & log2FoldChange > 0 ~ "UP",
      padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    ))
  
  return(labeled)
}



#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  plot <- ggplot(labeled_results, aes(x = pvalue)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    labs(
      title = "Distribution of Unadjusted P-values",
      x = "Unadjusted P-value",
      y = "Count"
    ) +
    theme_bw()
  
  return(plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  # Filter for significant genes based on padj threshold
  sig_results <- labeled_results %>%
    dplyr::filter(!is.na(padj) & padj < padj_threshold)
  
  # Build histogram of log2FoldChange values
  plot <- ggplot2::ggplot(sig_results, ggplot2::aes(x = log2FoldChange)) +
    ggplot2::geom_histogram(bins = 30,
                            fill = "steelblue",
                            color = "white",
                            alpha = 0.8) +
    ggplot2::labs(
      title = paste0("Distribution of Log2 Fold Changes"),
      subtitle = paste0("Genes significant at padj < ", padj_threshold),
      x = "Log2 Fold Change",
      y = "Count"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "dashed",
                        color = "red",
                        linewidth = 0.8)
  
  return(plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  
  # Get top genes ranked by ascending padj, removing NAs
  top_genes <- labeled_results %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = num_genes) %>%
    dplyr::pull(genes)  
  
  # Extract normalized counts from dds object for top genes
  norm_counts <- DESeq2::counts(dds_obj, normalized = TRUE)
  
  # Subset to top genes and pivot to long format for ggplot
  plot_data <- norm_counts[top_genes, ] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "sample",
      values_to = "normalized_counts"
    )
  
  # Build scatter plot
  plot <- ggplot2::ggplot(plot_data,
                          ggplot2::aes(x = gene,
                                       y = normalized_counts,
                                       color = sample)) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = paste0("Normalized Counts for Top ", num_genes, " Significant Genes"),
      subtitle = "Ranked by ascending adjusted p-value",
      x = "Gene",
      y = "Normalized Counts (log10 scale)",
      color = "Sample"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  return(plot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  
  plot_data <- labeled_results %>%
    dplyr::filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
    dplyr::mutate(neg_log10_padj = -log10(padj))
  
  plot <- ggplot2::ggplot(plot_data,
                          ggplot2::aes(x = log2FoldChange,
                                       y = neg_log10_padj,
                                       color = volc_plot_status)) +  # <-- not 'status'
    ggplot2::geom_point(size = 1.5, alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = c(
        "UP"   = "#E41A1C",
        "DOWN" = "#377EB8",
        "NS"   = "grey60"
      )
    ) +
    ggplot2::labs(
      title = "Volcano Plot of DESeq2 Results",
      x     = "Log2 Fold Change",
      y     = "-log10(Adjusted P-Value)",
      color = "Status"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title      = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    ) +
    ggplot2::geom_vline(xintercept = c(-1, 1),
                        linetype  = "dashed",
                        color     = "black",
                        linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(0.05),
                        linetype  = "dashed",
                        color     = "black",
                        linewidth = 0.5)
  
  return(plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  
  # Step 1: Load id2gene mapping
  id2gene <- readr::read_tsv(
    id2gene_path,
    col_names = c("ensembl_id", "gene_symbol"),
    show_col_types = FALSE
  )
  
  # Step 2: Filter out NA log2FoldChange before join
  filtered <- labeled_results %>%
    dplyr::filter(!is.na(log2FoldChange))
  
  # Step 3: Join gene symbols — do this as a standalone step
  joined <- dplyr::left_join(filtered, id2gene, by = c("genes" = "ensembl_id"))
  
  # Step 4: Verify gene_symbol exists after join
  stopifnot("gene_symbol" %in% colnames(joined))
  
  # Step 5: Filter, deduplicate, and rank
  ranked <- joined %>%
    dplyr::filter(!is.na(gene_symbol)) %>%
    dplyr::group_by(gene_symbol) %>%
    dplyr::slice_max(order_by = abs(log2FoldChange), n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(log2FoldChange))
  
  # Step 6: Build named vector
  rnk_list <- ranked$log2FoldChange
  names(rnk_list) <- ranked$gene_symbol
  
  return(rnk_list)
}
#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  # Load gene sets from GMT file
  pathways <- fgsea::gmtPathways(gmt_file_path)
  
  # Run fgsea with the ranked list and gene set size filters
  fgsea_res <- fgsea::fgsea(
    pathways  = pathways,
    stats     = rnk_list,
    minSize   = min_size,
    maxSize   = max_size
  )
  
  # Convert to tibble and arrange by ascending padj for readability
  fgsea_results <- fgsea_res %>%
    tibble::as_tibble() %>%
    dplyr::arrange(padj)
  
  return(fgsea_results)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths) {
  
  top_positive <- fgsea_results %>%
    dplyr::filter(!is.na(padj) & NES > 0) %>%
    dplyr::arrange(dplyr::desc(NES)) %>%
    dplyr::slice_head(n = num_paths)
  
  top_negative <- fgsea_results %>%
    dplyr::filter(!is.na(padj) & NES < 0) %>%
    dplyr::arrange(NES) %>%
    dplyr::slice_head(n = num_paths)
  
  plot_data <- dplyr::bind_rows(top_positive, top_negative) %>%
    dplyr::mutate(
      pathway = stringr::str_replace_all(pathway, "_", " "),
      pathway = stringr::str_trunc(pathway, width = 30, ellipsis = "..."),  # <-- truncate to 30 chars
      pathway = factor(pathway, levels = unique(pathway[order(NES)]))
    )
  
  plot <- ggplot2::ggplot(plot_data,
                          ggplot2::aes(x = NES,
                                       y = pathway,
                                       fill = NES > 0)) +
    ggplot2::geom_bar(stat = "identity", width = 0.7) +
    ggplot2::scale_fill_manual(
      values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
      labels = c("TRUE" = "Positive NES", "FALSE" = "Negative NES"),
      name   = "Direction"
    ) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "black") +
    ggplot2::labs(
      title = paste0("Top ", num_paths, " Pathways by Positive and Negative NES"),
      x     = "Normalized Enrichment Score (NES)",
      y     = "Pathway"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y     = ggplot2::element_text(size = 8, angle = 20, hjust = 1),  # <-- angled text
      axis.title      = ggplot2::element_text(face = "bold"),
      legend.position = "top"
    )
  
  return(plot)
}