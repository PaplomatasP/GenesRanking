#' Identify important genes in single-cell RNA-seq data using SCMarker
#'
#' This function uses the SCMarker package to identify important genes in single-cell RNA-seq data.
#'
#' @param data A matrix with rows representing genes and columns representing cells.
#' @param GeneSK An integer indicating the number of genes to use in the clustering.
#' @param CellSK An integer indicating the number of cells to use in the clustering.
#'
#' @return A list containing the SCMarker object and a vector of important genes.
#' @importFrom SCMarker ModalFilter GeneFilter getMarker
SCMarkerfun = function(data, GeneSK, CellSK) {


  obj <- as.matrix(data)


  res <-
    SCMarker::ModalFilter(
      data = obj,
      geneK = GeneSK,
      cellK = CellSK,
      cutoff = 3,
      width = 5
    )# default width = 1 for UMI data, width =2 for TPM data.
  res=SCMarker::getMarker(obj=res,k=300,n=30)
  genes = res[["geneSumm"]][["gene"]]
  Score = res[["geneSumm"]][["count"]]

  genesImp=res[["marker"]]

  FilterData <- data[ rownames(data) %in%  genesImp,]


  return(list(FilterData=FilterData,res=res[3:4] ,Important_Features=genesImp))

}


#' Select top self-expressed genes
#'
#' This function selects the top self-expressed genes based on the number of times a gene is its own nearest neighbor in a gene expression dataset.
#'
#' @param data A matrix or data frame where rows are samples and columns are genes.
#' @param n An integer indicating the number of top self-expressed genes to select.
#'
#' @return A list with two elements:
#' \describe{
#' \item{FilterData}{A data frame containing the top self-expressed genes.}
#' \item{Important_Features}{A character vector containing the names of the selected genes.}
#' }
#'
#' @importFrom SelfE SelfE
SelfEGenes = function(data, n) {


  obj <- as.data.frame(t(data)  )


  # Calling SelfE function
  GeneIDs <- SelfE::SelfE(obj, n)
  # Selecting genes from the data
  FilterData = obj[, GeneIDs]
  genes=colnames(FilterData)


  return(list(FilterData=FilterData,Important_Features=genes) )
}

#' DUBStepRfun
#'
#' This function filters cells using DUBStepR, which is an R package that enables identification of "drop-outs" in single-cell transcriptome analysis.
#'
#' @param data A data frame or matrix containing gene expression values.
#' @param K An integer indicating the number of nearest neighbors to use in determining cells that may be drop-outs.
#' @param Num.pcs An integer indicating the number of principal components to use for clustering.
#'
#' @return A list of the filtered gene expression data and a list of DUBStepR output results.
#'
#' @import DUBStepR
#' @import SingleCellExperiment
DUBStepRfun = function(data, K,Num.pcs) {
  obj = as.data.frame(data)


    sce <- SingleCellExperiment(assays = list(counts = obj))

    object_Seurat <-
      as.matrix(SummarizedExperiment::assay(sce, "counts"))


    Input <-
      Seurat::CreateSeuratObject(counts = object_Seurat)

    dubstepR.out <-
      DUBStepR::DUBStepR(
        input.data = Input@assays$RNA@data,
        min.cells = 0.05 * ncol(Input),
        optimise.features = T,
        k = K,
        num.pcs = Num.pcs,
        error = 0
      )
    Input@assays$RNA@var.features <-
      dubstepR.out$optimal.feature.genes
    genes = Input@assays$RNA@var.features



  FilterData = data[rownames(data) %in% genes, ]

  return(list(FilterData=FilterData,Important_Features=genes,dubstepR.out=dubstepR.out[["corr.info"]]) )
}




#' Perform scPNMF filtering on single-cell data
#'
#' @param data A numeric matrix of gene expression values
#' @param DM The distance metric to be used in the scPNMF algorithm (either "KL" for Kullback-Leibler divergence or "JS" for Jensen-Shannon divergence)
#' @param gN The number of genes to be retained in the filtered dataset
#'
#' @return A list with two elements:
#' \describe{
#'  \item{FilterData}{The filtered single-cell expression data}
#'  \item{Important_Features}{The selected genes with the highest relevance scores}
#' }
#'
#' @importFrom scPNMF PNMFfun basisSelect getInfoGene
scPNMFfun = function(data, DM,gN ) {


  obj=as.matrix(data)
  res_pnmf <- scPNMF::PNMFfun(X = as.matrix(data),
                              K = 3,
                              method=DM,
                              tol=1e-4,
                              maxIter=100,
                              verboseN = FALSE)




  W <- res_pnmf$Weight
  S <- res_pnmf$Score

  W_select <- scPNMF::basisSelect(
    W = W,
    S = S,
    X = obj,
    toTest = TRUE,
    toAnnotate = FALSE,
    mc.cores = 1
  )

  Score = as.data.frame(apply(W, 1, sum))
  W_select=as.data.frame(W_select)
  ig <-
    scPNMF::getInfoGene(
      W_select,
      M =  gN,
      by_basis = FALSE,
      return_trunW = TRUE,
      dim_use = NULL
    )


  FilterData = data[rownames(data) %in%  ig[["InfoGene"]], ]

  return(list(FilterData=FilterData,Important_Features=ig) )
}


#' Perform feature selection using M3Drop
#'
#' This function performs feature selection on a matrix of count data using M3Drop.
#' M3Drop is a method for detecting and filtering out cells with mitochondrial gene expression
#' patterns, which may be indicative of mitochondrial dysfunction or cell stress.
#'
#' @param data A matrix of count data.
#' @param Mt_Method The method to use for selecting mitochondrial genes. Must be one of "M3Drop", "percent", or "counts".
#' @param Mt_threshold The threshold to use for selecting mitochondrial genes. If Mt_Method is "M3Drop", this should be a numeric value between 0 and 1. If Mt_Method is "percent", this should be a percentage value between 0 and 100. If Mt_Method is "counts", this should be an integer value.
#' @return A list containing the filtered count data and a list of important features selected by M3Drop.
#'
#'
#' @importFrom M3Drop M3DropFeatureSelection
M3Dropfun = function(data,Mt_Method,Mt_threshold) {


      obj = as.matrix(data)

      DropsGenes <-
        M3Drop::M3DropFeatureSelection(
          obj,
          mt_method = Mt_Method,
          mt_threshold = Mt_threshold,
          suppress.plot = TRUE
        )



  FilterData = data[ rownames(data) %in%  rownames(DropsGenes),]


  return(list(FilterData=FilterData,Important_Features=rownames(DropsGenes),DropsGenes=DropsGenes) )
}

#' Differential expression analysis using the scDD method
#'
#' This function performs differential expression analysis on single-cell RNA sequencing data using the scDD method. It takes as input a count matrix, a list of prior parameters, and a vector of condition labels. The function returns a list of filtered data, important features, and results of the analysis.
#'
#' @param data A numeric matrix of raw count data.
#' @param prior_param A list of prior parameters used for Bayesian inference. Default values are alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01.
#' @param labels A vector of condition labels corresponding to the columns of the count matrix.
#' @param nGenes The number of top-ranked genes to return. Default is 150.
#'
#' @return A list with the following components:
#' \item{FilterData}{A filtered version of the input data, retaining only genes that were found to be differentially expressed.}
#' \item{Important_Features}{A vector of the names of the genes that were found to be differentially expressed.}
#' \item{Results}{A data frame with the results of the differential expression analysis.}
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom scDD scDD
scDDfun = function(data,prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01),labels,nGenes=150) {
  obj = as.matrix(data)

  condition <- labels

  sce <- SingleCellExperiment(assays= list(normcounts=obj, Counts=obj), colData=(condition))
  names(sce@colData@listData) <- "condition"

  ## ----main engine-----------------------------------------------------------
  scDatExSim <- scDD(sce, prior_param=prior_param, testZeroes=FALSE)

  ## ----main results----------------------------------------------------------
  RES <- scDD::results(scDatExSim)
  RES <- RES[order(RES$nonzero.pvalue.adj), ]

  Genes=RES$gene[1:nGenes]
  FilterData = data[ rownames(data) %in%  Genes,]

  return(list(FilterData=FilterData,Important_Features=Genes,Results=RES) )
}

#' Run SIMLR clustering and feature ranking on single-cell RNA-seq data
#'
#' This function runs the SIMLR algorithm on single-cell RNA-seq data to identify cell clusters
#' and rank the most informative genes. The function takes as input a matrix of gene expression
#' values, the number of clusters to identify, the ratio of cores to use for parallel processing,
#' and the number of top-ranked genes to return.
#'
#' @param data A matrix of gene expression values, where rows represent genes and columns represent cells.
#' @param ClusterNumber number of clusters to be estimated over your data.
#' that all available cores should be used.
#' @param nGenes The number of top-ranked genes to return. Default is 150.
#'
#' @return A list containing the filtered gene expression data, the list of top-ranked genes,
#' and the SIMLR feature ranking results.
#'
#' @importFrom SIMLR SIMLR SIMLR_Feature_Ranking
SIMLRFun=function(data,ClusterNumber,nGenes){

  RunSIMLR = SIMLR(X = data, c = ClusterNumber, cores.ratio = 0)
  Rank=SIMLR_Feature_Ranking(A=RunSIMLR[["S"]],X=data)
  Rank$aggR = rownames(data)[Rank$aggR]

  Genes=Rank$aggR[1:nGenes]
  FilterData = data[ rownames(data) %in%  Genes,]

  return(list(FilterData=FilterData,Important_Features=Genes,Results=Rank) )

}


#' Perform feature selection and identify important features using scmap
#'
#' This function performs feature selection and identifies important features based
#' on the output of the scmap package. The scmap package is used to find the most
#' informative genes for separating two groups of cells, and these genes are selected
#' as important features. The function returns a list containing the filtered data,
#' the list of important features, and the results of the scmap analysis.
#'
#' @param data A matrix or data frame containing the gene expression data for the single cells.
#' Rows represent genes, and columns represent cells.
#'
#' @param labels A vector containing the cell labels.
#' Rows represent cells, and columns represent metadata associated with the cells (e.g.,
#' sample ID, batch, cell type, etc.).
#'
#' @param nGenes An integer specifying the number of genes to select as important features.
#'
#' @return A list containing the following items:
#' \describe{
#'   \item{FilterData}{A matrix containing the filtered gene expression data. Rows represent genes,
#'   and columns represent cells. Only the nGenes genes selected as important features are included.}
#'   \item{Important_Features}{A character vector containing the names of the nGenes genes selected
#'   as important features.}
#'   \item{Results}{A data frame containing the output of the scmap analysis. Rows represent genes,
#'   and columns represent various statistics (e.g., scores, p-values, etc.) calculated by scmap.}
#' }
#'
#' @importFrom scmap selectFeatures
#' @importFrom SummarizedExperiment rowData
scmapfun=function(data,labels,nGenes){
  sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(data)), colData = labels)
  logcounts(sce) <- data
  SummarizedExperiment::rowData(sce)$feature_symbol <- rownames(sce)
  # remove features with duplicated names
  sce <- sce[!duplicated(rownames(sce)), ]


  sce <- selectFeatures(sce, suppress_plot = TRUE,n_features = nGenes)

  Res=as.data.frame(SummarizedExperiment::rowData(sce) )

  ResTRUE=subset(Res,Res$scmap_features==TRUE)
  ResTRUE=ResTRUE[order(ResTRUE$scmap_scores,decreasing = TRUE), ]

  Genes=rownames(ResTRUE)


  FilterData = data[ rownames(data) %in%  Genes,]

  return(list(FilterData=FilterData,Important_Features=Genes,Results=Res) )
}
