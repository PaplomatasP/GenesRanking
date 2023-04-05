#' Filter genes according to the specified method
#'
#' This function filters genes in scRNA-seq data according to the specified method.
#' The user can choose between two filtering methods: "High_Variable_Genes" and "Remove_Low_Variance".
#' If "High_Variable_Genes" is selected, the function uses the FindVariableFeatures function of Seurat to identify the top-ranked variable genes.
#' If "Remove_Low_Variance" is selected, the function removes genes with low variance using the nearZeroVar function of caret.
#'
#' @param data A numeric matrix of gene expression values with rows representing genes and columns representing cells.
#' @param filter_method A character string specifying the filtering method to be used. The valid options are "High_Variable_Genes" and "Remove_Low_Variance".
#' @param n_features An integer value that sets the number of top-ranked variable genes to keep. Only used when filter_method is "High_Variable_Genes".
#' @param unique_cut A numeric value that sets the threshold for low variance genes. Only used when filter_method is "Remove_Low_Variance".
#'
#' @return A filtered data frame with rows representing genes and columns representing cells.
#'
#' @examples
#'
#' data(ExampleDataset)
#' data_filtered <- filter_genes(ExampleDataset[1:100,],filter_method="High_Variable_Genes",n_features= 80)
#'
#' @importFrom Seurat CreateSeuratObject FindVariableFeatures VariableFeatures
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
#'
#' data(ExampleDataset)
#' Norm_data=normalize_data(ExpressionData=ExampleDataset[1:100,],normal_method="LogNormalize")
#'
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


