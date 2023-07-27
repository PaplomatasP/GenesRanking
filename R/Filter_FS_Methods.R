#' Filter Feature Selection Methods
#'
#' This function applies selected feature selection methods to single-cell RNA 
#' sequencing data. It takes as input a gene expression matrix, sample labels, 
#' the number of top genes to retain, and the name of the feature selection method. 
#' Optional parameters include log transformation and high variability gene filtering.
#'  The function returns a matrix of the top genes selected by the specified method.
#' Available feature selection methods include statistical tests 
#' (Wald Test, Beta Poisson GLM, Wilcoxon Test, Likelihood Ratio Test, Kruskal Wallis), 
#' correlation measures (Kendall, Spearman), information gain, and mRMRe.
#' This function is designed to facilitate the comparison of different feature 
#' selection methods and to assist in the identification of the most informative genes
#' for further analysis.
#'
#' @param data Data frame/matrix of gene expression values. Rows are genes and 
#'   columns are samples.
#' @param Labels Vector of sample labels. Length must equal the number of data columns.
#' @param n_genes_to_keep Integer. Number of top genes to retain after feature 
#'   selection.
#' @param method Character. Name of the feature selection method to apply. Valid 
#'   methods: "WaldTest", "BetaPoissonGLM", "WilcoxonTest", "LikelihoodRatioTest", 
#'   "KruskalWallis", "mRMRe", "QuasiLikelihoodFTest", "KendallCorrelation", 
#'   "SpearmanCorrelation", "InformationGain".
#' @param LogTransformation Logical. If TRUE, apply log transformation on the data 
#'   using Seurat. Default is TRUE.
#' @param HighVariableFIlter Logical. If TRUE, apply high variability gene filtering. 
#'   Default is TRUE.
#' @param n_features Integer. Number of genes with high variability to retain.
#' @return Data frame/matrix of top n_genes_to_keep genes selected by the specified 
#'   method.
#' @examples
#' # Example data and labels
#'
#' # Apply the Wald Test method and keep the top 10 genes
#' result <- Filter_FS_Methods(ExampleDataset[1:200, ], Labels,
#'     n_genes_to_keep = 50, method = "WaldTest",
#'     LogTransformation = TRUE, HighVariableFIlter = TRUE, n_features = 100)
#' @export
Filter_FS_Methods <- function(data,
                              Labels,
                              n_genes_to_keep,
                              method,
                              LogTransformation = TRUE,
                              HighVariableFIlter = TRUE,
                              n_features = 2000) {
    # Filter lowly-expressed genes
 min_genes <- 2000
min_cell_percentage <- 0.1
keep_going <- TRUE

while (keep_going) {
  # Filter lowly-expressed genes
  min_cells <- ceiling(min_cell_percentage * ncol(data)) # At least min_cell_percentage of cells
  gene_counts <- rowSums(data > 0)
  filtered_data <- data[gene_counts >= min_cells, ]
  
  # Check if we have at least 2000 genes
  if (nrow(filtered_data) >= min_genes) {
    keep_going <- FALSE
  } else {
    # If not, reduce min_cell_percentage to keep more genes
    min_cell_percentage <- min_cell_percentage - 0.01
  }
}

# Set data to the filtered data
data <- filtered_data

    if (HighVariableFIlter == TRUE) {
        seurat_obj <- Seurat::CreateSeuratObject(data)
        data_seurat <- Seurat::FindVariableFeatures(seurat_obj,
            selection.method = "vst",
            nfeatures = n_features
        )
        genes <- head(Seurat::VariableFeatures(data_seurat), n_features)
        data <- as.data.frame(data[rownames(data) %in% genes, ])
    }
    if (LogTransformation == TRUE) {
        scale_factor <- mean(colSums(data))
        data <-
            as.matrix(Seurat::LogNormalize(data, scale.factor = scale_factor))
    }

    result <- switch(method,
        "WaldTest" = Statistical_Filter(
            data,
            Pvalue_md = "Waldtest",
            Labels = Labels,
            n_genes_to_keep = n_genes_to_keep
        ),
        "BetaPoissonGLM" = Statistical_Filter(
            data,
            Pvalue_md = "BPglm",
            Labels = Labels,
            n_genes_to_keep = n_genes_to_keep
        ),
        "WilcoxonTest" = Statistical_Filter(
            data,
            Pvalue_md = "WilcoxonTest",
            Labels = Labels,
            n_genes_to_keep = n_genes_to_keep
        ),
        "LikelihoodRatioTest" = Statistical_Filter(
            data,
            Pvalue_md = "LRT",
            Labels = Labels,
            n_genes_to_keep = n_genes_to_keep
        ),
        "KruskalWallis" = KruskalWallisFun(data,
            Labels = Labels, n_genes_to_keep =
                n_genes_to_keep
        ),
        "mRMRe" = mRMReFun(data,
            Labels = Labels, n_genes_to_keep =
                n_genes_to_keep
        ),
        "QuasiLikelihoodFTest" = QLFTfun(data,
            Labels = Labels, n_genes_to_keep =
                n_genes_to_keep
        ),
        "KendallCorrelation" = KendallFun(data,
            Labels = Labels, n_genes_to_keep =
                n_genes_to_keep
        ),
        "SpearmanCorrelation" = CorrelationFun(data,
            Labels =
                Labels, n_genes_to_keep = n_genes_to_keep
        ),
        "InformationGain" = ig_fun(data,
            Labels = Labels, n_genes_to_keep =
                n_genes_to_keep
        ),
        stop(
            "Invalid method. Please choose a valid feature selection method."
        )
    )

    return(result)
}
