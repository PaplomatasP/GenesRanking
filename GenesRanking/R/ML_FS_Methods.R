#' ML Based Feature Selection Methods
#'
#' This function applies a specified machine learning method to filter and select important genes from single-cell RNA-seq data.
#'
#' @param data A matrix or data frame containing the input data (features and samples).
#' @param Labels A vector containing the true class labels of the samples.
#' @param n_genes_to_keep An integer representing the number of top-ranked genes to retain based on the selected method.
#' @param method A character string specifying the machine learning method to use for feature selection. Supported methods:
#'  "PAM", "gcvEarth", "PLS", "RF", "CART", "C5.0", "GBM", "xgbTree", "Catboost", and "ShapleyValues".
#' @param LogTransformation Logical, if TRUE log transformation is applied on the data using the Seurat package. Default is TRUE.
#' @param HighVariableFIlter Logical, if TRUE high variability gene filtering is applied. Default is TRUE.
#' @param n_features Integer, the number of genes with high variability to keep.
#'
#' @return A list containing the filtered data, important features, and any additional information relevant to the selected method.
#' @examples
#' # Example usage:
#' data(ExampleDataset)
#' data(Labels)
#' # Apply the knn method and keep the top 10 genes
#' result <- ML_FS_Methods(ExampleDataset[1:500,], Labels, n_genes_to_keep=100, method="C5.0",
#' LogTransformation=TRUE, HighVariableFIlter=TRUE, n_features=200)
#' @export
ML_FS_Methods  <- function(data, Labels, n_genes_to_keep, method,LogTransformation=TRUE,
                           HighVariableFIlter=TRUE,n_features=2000) {
  min_cell_percentage<-0.1
  # Filter lowly-expressed genes
  min_cells <- ceiling(min_cell_percentage * ncol(data))  # At least min_cell_percentage of cells
  gene_counts <- rowSums(data > 0)
  data <- data[gene_counts >= min_cells,]

  if (HighVariableFIlter==TRUE){
    seurat_obj <- Seurat::CreateSeuratObject(data)
    data_seurat <- Seurat::FindVariableFeatures(
      seurat_obj,
      selection.method = "vst",
      nfeatures = n_features
    )
    genes <- head(Seurat::VariableFeatures(data_seurat), n_features)
    data <- as.data.frame(data[rownames(data) %in% genes, ])
  }
  if (LogTransformation==TRUE){
    scale_factor <- mean(colSums(data))
    data <- as.matrix(
      Seurat::LogNormalize(data, scale.factor <- scale_factor))
  }

  if (method == "PAM") {
    result <- ML_filter(data, MLmethod="pam", Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else if (method == "gcvEarth") {
    result <-  ML_filter(data, MLmethod="gcvEarth", Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else if (method == "PLS") {
    result <- ML_filter(data, MLmethod="pls", Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else if (method == "RF") {
    result <-ML_filter(data, MLmethod="rf", Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else if (method == "CART") {
    result <- ML_filter(data, MLmethod="rpart", Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else if (method == "C5.0") {
    result <-ML_filter(data, MLmethod="C5.0", Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else if (method == "GBM") {
    result <- ML_filter(data, MLmethod="gbm", Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else if (method == "xgbTree") {
    result <- ML_filter(data, MLmethod="xgbTree", Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else if (method == "Catboost") {
    result <- CatBoost_feature_selection(data,Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else if (method == "ShapleyValues") {
    result <- ShapleyValuesFun(data, Labels=Labels, n_genes_to_keep=n_genes_to_keep)
  } else {
    stop("Invalid method. Please choose a valid feature selection method.")
  }
  return(result)
}

