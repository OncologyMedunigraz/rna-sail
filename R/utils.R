#' RNASAIL: Comprehensive RNA-seq Analysis Pipeline
#'
#' This package provides functions for differential expression analysis,
#' pathway analysis, WGCNA co-expression networks, and visualization of RNA-seq data.
#'
#' @docType package
#' @name rna-sail

# Package startup message
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("RNASAIL loaded successfully!")
  packageStartupMessage("Use ?RNASAIL to get started.")

  # Set up conflicts preferences
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflicts_prefer(stats::var)
    conflicted::conflicts_prefer(dplyr::select)
  }
}

#' Install Missing Packages
#'
#' Helper function to install missing CRAN packages.
#'
#' @param packages Character vector of package names to check and install
#' @return None (invisible)
#' @export
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
  invisible()
}

#' Install Missing Bioconductor Packages
#'
#' Helper function to install missing Bioconductor packages.
#'
#' @param packages Character vector of Bioconductor package names
#' @return None (invisible)
#' @export
bioc_install_if_missing <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
    library(pkg, character.only = TRUE)
  }
  invisible()
}

#' Check Column Validity for PCA
#'
#' Helper function to check if a column has finite values and non-zero variance.
#'
#' @param x Numeric vector to check
#' @return Logical indicating if column is valid
#' @keywords internal
is_good_column <- function(x) {
  all(is.finite(x)) && var(x) > 0
}

#' Create Output Directory
#'
#' Creates output directory if it doesn't exist.
#'
#' @param dir_path Path to directory to create
#' @return Character path to directory
#' @export
create_output_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created directory: ", dir_path)
  }
  return(dir_path)
}

#' Get Default Required Packages
#'
#' Returns the default list of required packages for the pipeline.
#'
#' @return Character vector of package names
#' @export
get_required_packages <- function() {
  c(
    "rtracklayer", "limma", "edgeR", "fgsea", "msigdbr",
    "dplyr", "ggplot2", "ggrepel", "RColorBrewer", "colorspace",
    "ComplexHeatmap", "circlize", "data.table", "stringr",
    "conflicted", "scales", "tidyverse", "matrixStats"
  )
}

#' Get Default Bioconductor Packages
#'
#' Returns the default list of required Bioconductor packages.
#'
#' @return Character vector of Bioconductor package names
#' @export
get_required_bioc_packages <- function() {
  c(
    "limma", "edgeR", "msigdbr", "DOSE", "clusterProfiler",
    "WGCNA", "viper", "dorothea", "rtracklayer", "org.Mm.eg.db",
    "org.Hs.eg.db", "biomaRt"
  )
}


#' Check if a series of Packages is installed
#'
#' @param required_pkgs vector of names of the required packages
#' @return Prints uninstalled packages
#' @export
check_packages <- function(required_pkgs) {
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required but not installed")
    }
  }
}


#' P-value Calculator for the fn plot_gene_expression()
#'
#' @param data a dataframe with a column containing the gene expr values and more
#' @param facet.col name of the column according to which group
#' @param y.col name of the column with gene data
#' @param x.col name of the column of inside grouping
#' @param stat_test the kind of stat_test to perform
#' @param p_correction the type of p_correction to implement (default: "none")
#'
#' @return df with the p-vals for the various comparisons
#' @export
PvalCalc <- function(data, facet.col, y.col, x.col, stat_test, p_correction = "none") {
  
  CheckPackages(c("dplyr", "ggpubr"))
  
  comps <- combn(levels(factor(data[[x.col]])), m = 2, simplify = FALSE)
  
  res <- ggpubr::compare_means(
    formula    = as.formula(paste(y.col, "~", x.col)),
    data       = data,
    method     = stat_test,
    comparisons = comps,
    p.adjust.method = p_correction,
    group.by   = facet.col
  )
  
  # Calculate y positions row-wise using mapply
  # (avoids coercion to character matrix)
  res$y.position <- mapply(function(g1, g2, facet_val) {
    max(data[data[[x.col]] %in% c(g1, g2) & data[[facet.col]] == facet_val, y.col])
  }, res$group1, res$group2, res[[facet.col]])
  
  # Offset y positions per facet to avoid label overlap
  res <- res %>%
    dplyr::group_by(.data[[facet.col]]) %>%
    dplyr::mutate(y.position = y.position + seq(0, by = 0.1, length.out = dplyr::n())) %>%
    dplyr::ungroup()
  
  return(res)
}



get_palette <- function(n_conditions, val_conditions) {
  check_packages("grDevices", "RColorBrewer")
  
  if (n_conditions <= 9) {
    color_mapping <- setNames(
      RColorBrewer::brewer.pal(n_conditions, "Set1"),
      val_conditions
    )
  } else {
    # !!! Must check grDevices is installed !!!
    pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    color_mapping <- setNames(
      pal(n_conditions),
      val_conditions
    )
  }

  return(color_mapping)

}
