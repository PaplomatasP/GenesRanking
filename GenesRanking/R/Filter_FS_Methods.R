#' Filter Feature Selection Methods
#'
#' This function applies various feature selection methods for single-cell RNA sequencing data.
#'
#' @param data A data frame or matrix containing gene expression values. Rows represent genes and columns represent samples.
#' @param Labels A vector containing the labels of the samples. Must be the same length as the number of columns in the data.
#' @param n_genes_to_keep Integer, the number of top genes to keep after applying the selected feature selection method.
#' @param method Character, the name of the feature selection method to be applied. Valid methods include:
#'   "WaldTest", "BetaPoissonGLM", "WilcoxonTest", "LikelihoodRatioTest", "KruskalWallis", "mRMRe",
#'   "QuasiLikelihoodFTest", "KendallCorrelation", "SpearmanCorrelation", and "InformationGain".
#' @param LogTransformation Logical, if TRUE log transformation is applied on the data using the Seurat package. Default is TRUE.
#' @param HighVariableFIlter Logical, if TRUE high variability gene filtering is applied. Default is TRUE.
#' @param n_features Integer, the number of genes with high variability to keep.
#' @return A data frame or matrix containing the top n_genes_to_keep genes selected by the specified method.
#'
#' @examples
#' # Example data and labels
#' data(ExampleDataset)
#' data(Labels)
#'
#' # Apply the Wald Test method and keep the top 10 genes
#' result <- Filter_FS_Methods(ExampleDataset[1:500, ], Labels,
#'     n_genes_to_keep = 100, method = "WaldTest",
#'     LogTransformation = TRUE, HighVariableFIlter = TRUE, n_features = 200
#' )
#' @export
Filter_FS_Methods <- function(data,
                              Labels,
                              n_genes_to_keep,
                              method,
                              LogTransformation = TRUE,
                              HighVariableFIlter = TRUE,
                              n_features = 2000) {
    min_cell_percentage <- 0.1
    # Filter lowly-expressed genes
    min_cells <-
        ceiling(min_cell_percentage * ncol(data)) # At least min_cell_percentage of cells
    gene_counts <- rowSums(data > 0)
    data <- data[gene_counts >= min_cells, ]
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
