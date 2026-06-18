#' Run LINCS Connectivity Analysis
#'
#' Performs connectivity mapping using LINCS L1000 data to identify compounds
#' that mimic or reverse the observed gene expression signature.
#'
#' @param de_results Differential expression results from limma
#' @param species Species ("mouse" or "human", default: "mouse")
#' @param n_up Number of top upregulated genes to use (default: 150)
#' @param n_down Number of top downregulated genes to use (default: 150)
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#'
#' @return List containing LINCS connectivity results
#'
#' @importFrom limma topTable
#' @importFrom AnnotationDbi mapIds
#' @importFrom utils head
#' @importFrom stats na.omit
#'
#' @export
run_lincs_analysis <- function(de_results, species = "mouse", n_up = 150, n_down = 150,
                               output_dir, experiment_name) {

  # Check required packages
  required_pkgs <- c("signatureSearch", "homologene")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed. ",
           "Install with: BiocManager::install('", pkg, "')")
    }
  }

  create_output_dir(output_dir)
  message("Starting LINCS connectivity analysis...")

  # Prepare differential expression data
  if (is.data.frame(de_results)) {
    deg_data <- de_results
  } else {
    # Extract from limma efit object
    deg_data <- topTable(de_results, number = Inf, adjust.method = "BH")
  }
  message("ok1")

  # Map mouse genes to human if needed
  print(deg_data)
  if (species == "mouse") {
    gene_mapping <- map_mouse_to_human_genes(deg_data, output_dir, experiment_name)
    up_genes_human <- gene_mapping$up_genes_human
    down_genes_human <- gene_mapping$down_genes_human

  } else {
    # rownames are SYMBOL_ENSG...
    symbols <- sub("_.*", "", rownames(deg_data))  # keep part BEFORE underscore

    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      stop("Package 'org.Hs.eg.db' is required for human gene mapping")
    }

    entrez_ids <- mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys     = symbols,
      keytype  = "SYMBOL",
      column   = "ENTREZID",
      multiVals = "first"
    )

    deg_data$entrez <- entrez_ids

    # Select top up and down genes (you can tweak thresholds if you want)
    up_genes_human   <- head(na.omit(deg_data[deg_data$logFC > 0,  "entrez"]), n_up)
    down_genes_human <- head(na.omit(deg_data[deg_data$logFC < 0,  "entrez"]), n_down)
  }

  message("ok2")
  # Check if we have enough genes
  if (length(up_genes_human) < 10 || length(down_genes_human) < 10) {
    warning("Insufficient genes for LINCS analysis. Need at least 10 up and 10 down genes.")
    return(NULL)
  }

  message("Using ", length(up_genes_human), " upregulated and ",
          length(down_genes_human), " downregulated human genes")

  # Create query signature
  query_signature <- create_lincs_query_signature(up_genes_human, down_genes_human)

  # Run LINCS connectivity analysis
  lincs_results <- perform_lincs_gess(query_signature, output_dir, experiment_name)

  # Create visualizations
  if (!is.null(lincs_results)) {
    create_lincs_visualizations(lincs_results, output_dir, experiment_name)
  }

  message("LINCS connectivity analysis completed!")

  return(lincs_results)
}




#' Map Mouse to Human Genes
#'
#' Maps mouse gene IDs to human orthologs using homologene.
#'
#' @param deg_data Differential expression data
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#'
#' @return List with human gene IDs
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom utils head write.table
#' @importFrom stats na.omit
#' @importFrom homologene homologene
#'
#' @keywords internal
map_mouse_to_human_genes <- function(deg_data, output_dir, experiment_name) {

  if (!requireNamespace("homologene", quietly = TRUE) ||
      !requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop("Packages 'homologene' and 'org.Mm.eg.db' are required for mouse-to-human mapping")
  }

  message("Mapping mouse genes to human orthologs...")

  # Extract gene symbols (assuming format is SYMBOL_ENSEMBLID)
  symbols <- sub(".*_", "", rownames(deg_data))

  # Map mouse symbols to Entrez IDs
  mouse_entrez <- mapIds(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = symbols,
    keytype = "ENSEMBL",
    column = "ENTREZID",
    multiVals = "first"
  )

  deg_data$mouse_entrez <- mouse_entrez

  # Select top up and down genes
  up_genes_mouse <- head(na.omit(deg_data[deg_data$logFC > 1, "mouse_entrez"]), 200)
  down_genes_mouse <- head(na.omit(deg_data[deg_data$logFC < -1, "mouse_entrez"]), 200)

  # Map to human using homologene
  mouse_genes_all <- unique(c(up_genes_mouse, down_genes_mouse))
  mouse_genes_numeric <- as.numeric(mouse_genes_all[!is.na(mouse_genes_all)])

  # Perform homolog mapping
  homolog_map <- homologene(mouse_genes_numeric, inTax = 10090, outTax = 9606)

  if (nrow(homolog_map) == 0) {
    stop("No homologous genes found between mouse and human")
  }

  # Map up and down genes separately
  up_human_map <- homolog_map$`9606_ID`[homolog_map$`10090_ID` %in% as.numeric(up_genes_mouse)]
  down_human_map <- homolog_map$`9606_ID`[homolog_map$`10090_ID` %in% as.numeric(down_genes_mouse)]

  # Convert to character and remove duplicates
  up_genes_human <- unique(as.character(up_human_map))
  down_genes_human <- unique(as.character(down_human_map))

  # Save mapping results
  mapping_file <- file.path(output_dir, paste0(experiment_name, "_mouse_human_gene_mapping.tsv"))
  mapping_df <- data.frame(
    mouse_entrez = homolog_map$`10090_ID`,
    human_entrez = homolog_map$`9606_ID`
  )
  write.table(mapping_df, file = mapping_file, sep = "\t", quote = FALSE, row.names = FALSE)

  message("Mapped ", length(mouse_genes_numeric), " mouse genes to ", nrow(homolog_map), " human orthologs")
  message("Final signature: ", length(up_genes_human), " up, ", length(down_genes_human), " down")

  return(list(
    up_genes_human = up_genes_human,
    down_genes_human = down_genes_human,
    mapping = homolog_map
  ))
}



#' Create LINCS Query Signature
#'
#' Creates a query signature object for LINCS analysis.
#'
#' @param up_genes Character vector of upregulated gene IDs
#' @param down_genes Character vector of downregulated gene IDs
#'
#' @return signatureSearch qSig object
#'
#' @importFrom signatureSearch qSig
#'
#' @keywords internal
create_lincs_query_signature <- function(up_genes, down_genes) {

  if (!requireNamespace("signatureSearch", quietly = TRUE)) {
    stop("Package 'signatureSearch' is required")
  }

  message("Creating LINCS query signature...")

  # Create query signature
  query_sig <- qSig(
    query = list(upset = up_genes, downset = down_genes),
    gess_method = "LINCS",
    refdb = "lincs"
  )

  return(query_sig)
}





#' Perform LINCS GESS Analysis
#'
#' Runs the actual LINCS gene expression signature search.
#'
#' @param query_signature Query signature object
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#'
#' @return LINCS results object
#'
#' @importFrom signatureSearch gess_lincs result
#' @importFrom utils write.table
#'
#' @keywords internal
perform_lincs_gess <- function(query_signature, output_dir, experiment_name) {

  if (!requireNamespace("signatureSearch", quietly = TRUE)) {
    stop("Package 'signatureSearch' is required")
  }

  # Ensure the data object exists
  if (!"clue_moa_list" %in% ls(envir = .GlobalEnv)) {

    data("clue_moa_list", package = "signatureSearch", envir = environment())
  }

  message("Running LINCS connectivity search... (this may take several minutes)")

  tryCatch({
    lincs_results <- gess_lincs(
      qSig    = query_signature,
      sortby  = "NCS",
      tau     = TRUE
    )

    results_table <- result(lincs_results)

    results_file <- file.path(output_dir, paste0(experiment_name, "_LINCS_results.tsv"))
    write.table(results_table, file = results_file, sep = "\t",
                quote = FALSE, row.names = FALSE)

    message("LINCS analysis completed. Found ", nrow(results_table), " connections")
    message("Results saved to: ", results_file)

    return(list(
      lincs_object  = lincs_results,
      results_table = results_table
    ))

  }, error = function(e) {
    message("LINCS analysis failed: ", e$message)
    message("This is likely a missing data object from signatureSearchData, ",
            "or (less likely) a connectivity issue.")
    return(NULL)
  })
}




#' Create LINCS Visualizations
#'
#' Creates various visualizations for LINCS connectivity results.
#'
#' @param lincs_results LINCS results list
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#'
#' @return None (creates plots)
#'
#' @keywords internal
create_lincs_visualizations <- function(lincs_results, output_dir, experiment_name) {

  results_table <- lincs_results$results_table
  lincs_object <- lincs_results$lincs_object

  # Create connectivity score plot
  create_lincs_connectivity_plot(results_table, output_dir, experiment_name)

  # Create drug class analysis
  create_lincs_drug_class_plot(results_table, output_dir, experiment_name)

  # Create top compounds visualization (if signatureSearch supports it)
  if (!is.null(lincs_object)) {
    create_lincs_heatmap(lincs_object, results_table, output_dir, experiment_name)
  }
}




#' Create LINCS Connectivity Score Plot
#'
#' Creates a bar plot showing connectivity scores for top compounds.
#'
#' @param results_table LINCS results table
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @param n_compounds Number of compounds to show (default: 30)
#'
#' @return ggplot object
#'
#' @import ggplot2
#' @importFrom stats reorder
#'
#' @keywords internal
create_lincs_connectivity_plot <- function(results_table, output_dir, experiment_name, n_compounds = 30) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Skipping LINCS connectivity plot - ggplot2 not available")
    return(NULL)
  }

  message("Creating LINCS connectivity plot...")

  top_hits <- head(results_table, n_compounds)

  p <- ggplot(top_hits, aes(
    x = reorder(.data$pert, .data$NCS),
    y = .data$NCS,
    fill = .data$NCS
  )) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(
      title = paste("Top", n_compounds, "LINCS Connectivity Scores"),
      x = "Compound",
      y = "Normalized Connectivity Score (NCS)"
    ) +
    theme_minimal()

  output_file <- file.path(output_dir, paste0(experiment_name, "_LINCS_connectivity_plot.pdf"))
  ggsave(output_file, p, width = 8, height = 6)
  message("Saved connectivity plot to: ", output_file)

  return(p)
}




#' Create LINCS Drug Class Plot
#'
#' Summarizes connectivity scores by drug class if available.
#'
#' @param results_table LINCS results table
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#'
#' @return ggplot object
#'
#' @import ggplot2
#' @importFrom stats aggregate reorder
#'
#' @keywords internal
create_lincs_drug_class_plot <- function(results_table, output_dir, experiment_name) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("Skipping drug class plot - ggplot2 not available")
    return(NULL)
  }

  if (!"pert_type" %in% colnames(results_table)) {
    message("No drug class (pert_type) column found in results. Skipping drug class plot.")
    return(NULL)
  }

  message("Creating drug class summary plot...")

  class_summary <- aggregate(NCS ~ pert_type, data = results_table, FUN = mean, na.rm = TRUE)
  class_summary <- class_summary[order(-class_summary$NCS), ]

  p <- ggplot(class_summary, aes(
    x = reorder(.data$pert_type, .data$NCS),
    y = .data$NCS,
    fill = .data$NCS
  )) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(
      title = "Average LINCS Connectivity Score by Drug Class",
      x = "Drug Class",
      y = "Mean NCS"
    ) +
    theme_minimal()

  output_file <- file.path(output_dir, paste0(experiment_name, "_LINCS_drug_class_plot.pdf"))
  ggsave(output_file, p, width = 8, height = 6)
  message("Saved drug class plot to: ", output_file)

  return(p)
}




#' Create LINCS Heatmap
#'
#' Creates a heatmap of top compounds' connectivity signatures.
#'
#' @param lincs_object LINCS analysis object
#' @param results_table Results table
#' @param output_dir Output directory
#' @param experiment_name Experiment name
#' @param n_compounds Number of compounds to show (default: 20)
#' @param n_genes Number of genes to show (default: 100)
#'
#' @return pheatmap object
#'
#' @importFrom utils head
#' @importFrom stats var
#' @importFrom signatureSearch getSig
#' @importFrom methods slot
#' @importFrom pheatmap pheatmap
#' @import grDevices
#'
#' @keywords internal
create_lincs_heatmap <- function(lincs_object, results_table, output_dir, experiment_name,
                                 n_compounds = 20, n_genes = 100) {

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message("Skipping heatmap - pheatmap not available")
    return(NULL)
  }

  if (!all(c("pert", "cell") %in% colnames(results_table))) {
    message("Results table must contain 'pert' and 'cell' columns. Skipping heatmap.")
    return(NULL)
  }

  message("Creating heatmap of top LINCS hits...")

  top_hits <- head(results_table, n_compounds)
  top_hits <- unique(top_hits[, c("pert", "cell")])


  refdb <- slot(lincs_object, "refdb")

  sig_mat <- try(
    getSig(
      cmp = top_hits$pert,
      cell = top_hits$cell,
      refdb = refdb
    ),
    silent = TRUE
  )

  if (inherits(sig_mat, "try-error") || is.null(sig_mat) || ncol(sig_mat) == 0) {
    message("Could not retrieve LINCS signatures for top hits. Skipping heatmap.")
    return(NULL)
  }


  sig_mat <- as.matrix(sig_mat)

  if (ncol(sig_mat) < 2 || nrow(sig_mat) < 2) {
    message("Signature matrix too small for heatmap. Skipping.")
    return(NULL)
  }


  sample_labels <- paste(top_hits$pert, top_hits$cell, sep = " | ")
  if (length(sample_labels) == ncol(sig_mat)) {
    colnames(sig_mat) <- make.unique(sample_labels)
  }

  gene_var <- apply(sig_mat, 1, var, na.rm = TRUE)
  gene_var[is.na(gene_var)] <- 0

  keep_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(n_genes, nrow(sig_mat)))]
  sig_sub <- sig_mat[keep_genes, , drop = FALSE]

  plot_mat <- t(sig_sub)
  output_file <- file.path(output_dir, paste0(experiment_name, "_LINCS_heatmap.pdf"))
  pdf(output_file, width = 10, height = 7)
  pheatmap(plot_mat,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           main = "Top LINCS Hits: Reference Signature Heatmap")
  dev.off()

  message("Saved heatmap to: ", output_file)
  invisible(plot_mat)
}
