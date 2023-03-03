#' Apply a statistical filter to gene expression data
#'
#' This function performs a statistical filtering method on gene expression data. It filters the genes based on their statistical significance calculated by various methods such as Seurat, DESeq2, MAST, or BPSC.
#'
#' @param data A matrix or data.frame object representing gene expression data. Rows represent genes, and columns represent samples.
#' @param filter_method A character string indicating the method to use for filtering. Can be one of "High_Variable_Genes" (to use the Seurat FindVariableFeatures method) or "Remove_Low_Variance" (to use the zero variance filter).
#' @param n_features An integer indicating the number of genes to keep after filtering using the Seurat FindVariableFeatures method. Default is 500.
#' @param unique_cut A numeric cutoff threshold for the zero variance filter. Features with a variance below this threshold will be filtered out.
#'
#' @return A filtered dataset with only the genes that pass the statistical filter.
#' \item{filtered_data}{A matrix or data.frame object representing the filtered gene expression data.}
#' \item{important_genes}{A character vector of gene names that pass the statistical filter.}
#' \item{statistical_analysis}{A data.frame object containing the p-values or FDRs for all genes, as calculated by the selected statistical method.}
#'
#' @examples
#' \donttest{
#' data(ExampleDataset)
#' data_filtered <- filter_genes(ExampleDataset,filter_method="High_Variable_Genes",n_features= 100)
#'}
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom BPSC BPglm
#' @importFrom MAST FromMatrix zlm lrTest
#' @importFrom utils head
#' @import Seurat
#' @import dplyr
#' @import tidyr
#' @importFrom caret nearZeroVar
#' @export
filter_genes = function(data, filter_method, n_features, unique_cut) {
  # Filter genes according to the specified method.
  if (filter_method == "High_Variable_Genes") {
    # Filter according to the FindVariableFeatures function of Seurat.
    seurat_obj <- Seurat::CreateSeuratObject(data)
    data_seurat <- Seurat::FindVariableFeatures(seurat_obj,
                                                selection.method = "vst",
                                                nfeatures = n_features)
    genes <- head(Seurat::VariableFeatures(data_seurat), n_features)
    filtered_data <- as.data.frame(data[rownames(data) %in% genes, ])
  } else if (filter_method == "Remove_Low_Variance") {
    # Filter according to the zero variance filter using a cutoff threshold.
    filtered_data <- data[-caret::nearZeroVar(data, uniqueCut = unique_cut), ]
  }
  return(filtered_data)
}

#' Normalize a matrix of gene expression data using different methods
#'
#' This function normalizes a matrix of gene expression data using different methods, including TMM, RLE, CPM, TPM, and LogNormalize.
#'
#' @param ExpressionData A Table of gene expression data to be normalized. Rows represent genes, and columns represent samples.
#' @param normal_method A character string specifying the normalization method. Possible values are "TMM", "RLE", "CPM", "TPM", and "LogNormalize".
#'
#' @return A matrix of normalized gene expression data.
#'
#' @examples
#' \donttest{
#' data(ExampleDataset)
#' Norm_data=normalize_data(ExpressionData=ExampleDataset,normal_method="LogNormalize")
#'}
#' @importFrom edgeR calcNormFactors
#' @importFrom stats median model.matrix na.omit p.adjust
#' @importFrom scater calculateTPM
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Seurat LogNormalize
#' @importFrom utils head
#' @export
normalize_data = function(ExpressionData, normal_method) {

  if(normal_method == "TMM") {
    norm_factor <-edgeR::calcNormFactors(ExpressionData,method=normal_method)
    normalized_data <- t(t(ExpressionData) / norm_factor)
  } else if(normal_method == "RLE") {
    geomeans <- exp(rowMeans(log(ExpressionData)))
    SF <- function(cnts) {
      norm_factor <- stats::median((cnts / geomeans)[is.finite(geomeans) & geomeans > 0])
    }
    norm_factor <- apply(ExpressionData, 2, SF)
    normalized_data <- t(t(ExpressionData) / norm_factor)
  } else if(normal_method == "CPM") {
    norm_factor <- colSums(ExpressionData)
    normalized_data <- t(t(ExpressionData) / norm_factor) * 1e6
  } else if(normal_method == "TPM") {
    Matrix_data <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = ExpressionData))
    normalized_data <- scater::calculateTPM(Matrix_data, exprs_values = "counts")
  } else if(normal_method == "LogNormalize") {
    scale_factor <- mean(colSums(ExpressionData))
    normalized_data <- as.matrix(Seurat::LogNormalize(ExpressionData,
                                            scale.factor = scale_factor) )
  }
  return(normalized_data)
}


