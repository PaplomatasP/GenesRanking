#' Variable_Filter function for feature selection in scRNA-seq data
#'
#' This function allows the user to apply various feature selection methods to single-cell RNA sequencing (scRNA-seq) data.
#' The function takes in a data matrix, a vector of cell labels, a feature selection method, and any additional arguments required by the method.
#' The resulting filtered data and important features are returned as a list.
#'
#' @param data A numeric matrix of gene expression values with rows representing genes and columns representing cells.
#' @param method A character string specifying the feature selection method to be used. The valid options are "SCMarker", "SelfE", "DUBStepR", "scPNMF", or "M3Drop".
#' @param ... Additional arguments that are required for the specific feature selection method chosen.
#' @return A list with two elements: \code{FilterData} and \code{Important_Features}.
#' \code{FilterData} is a numeric matrix of gene expression values with rows representing genes and columns representing cells, after the feature selection process.
#' \code{Important_Features} is a vector of important genes selected by the feature selection method.
#' @importFrom Seurat CreateSeuratObject
#' @importFrom DUBStepR DUBStepR
#' @importFrom M3Drop M3DropFeatureSelection
#' @importFrom scPNMF PNMFfun basisSelect getInfoGene
#'
#' @examples
#' \donttest{
#' data(ExampleDataset)
#' result <- Variable_Filter(ExampleDataset,method = "DUBStepR",K=10,Num.pcs=10)
#'}
#' @export
Variable_Filter <- function(data,method, ...) {

  if (method == "SCMarker") {
    result <- SCMarkerfun(data, GeneSK=10, CellSK=10)

  } else if (method == "SelfE") {
    result <- SelfEGenes(data,n=150)

  } else if (method == "DUBStepR") {
    result <- DUBStepRfun(data,K=10,Num.pcs=20)

  } else if (method == "scPNMF") {
    result <- scPNMFfun(data,DM="EucDist",gN=100 )

  } else if (method == "M3Drop") {
    result <- M3Dropfun(data,Mt_threshold=0.05,Mt_Method="fdr")

  }else if (method == "scDD") {
    result <- scDDfun(data,prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01),Labels,nGenes=150)

  }else if (method == "SIMLR") {
    result <- scDDfun(data,ClusterNumber=3,nGenes=150)

  }else if (method == "scmap") {
    result <- scmapfun(data,Labels=Labels,nGenes=150)
  }else {
    stop("Invalid method selected.")
  }


  return(result)
}


