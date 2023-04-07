#' Apply a statistical filter to gene expression data
#'
#' This function performs a statistical filtering method on gene expression data.
#' It filters the genes based on their statistical significance calculated by
#' various methods.
#'
#' @param data A matrix or data.frame object representing gene expression data.
#'   Rows represent genes, and columns represent samples.
#' @param Pvalue_md A character string indicating the statistical method to use
#'   for filtering. Can be one of "Waldtest", "BPglm", "Wilcoxon-Test" or "LRT".
#' @param Labels A character or factor vector specifying the label for each cell
#'   in \code{data}.
#' @param n_genes_to_keep An integer indicating the number of genes to keep after
#'   filtering.
#'
#' @return A list with the following components:
#' \item{FilterData}{A matrix or data.frame object representing the filtered
#'   gene expression data.}
#' \item{Important_Features}{A character vector of gene names that pass the
#'   statistical filter.}
#' \item{Statistical_Analysis}{A data.frame object containing the p-values and
#'   FDRs for all genes, as calculated by the selected statistical method.}
#'
#' @examples
#'
#' data(ExampleDataset)
#' Statistical_Filter(head(ExampleDataset, 10), Pvalue_md<-"Waldtest",
#' Labels<-Labels,n_genes_to_keep<-100)
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom BPSC BPglm
#' @importFrom MAST FromMatrix zlm lrTest
#' @importFrom stats model.matrix na.omit p.adjust
#' @importFrom utils head
#' @import Seurat
#'
#' @export
Statistical_Filter <- function(data, Pvalue_md, Labels, n_genes_to_keep) {
  options(scipen = 999)
  threshold <- 1
  data1 <- data_fitting(data, Labels)
  PvalueData <- PvalueCalc(data <- data1,
                          Labels <- Labels, Pvaluemethod = Pvalue_md)
  PvalueData  <- as.data.frame(PvalueData)
  if (Pvalue_md == "Waldtest" |
      Pvalue_md == "BPglm" | Pvalue_md == "WilcoxonTest") {
    rownames(PvalueData) <- PvalueData[, 1]
    count <- 0
    PvalueTreshold <- data.frame(Genes = character(),
                                Fdr = numeric(),
                                stringsAsFactors = FALSE)
    for (i in seq_len(nrow(PvalueData))) {
      if (PvalueData[, 3][i] <= threshold) {
        count <- count + 1
        PvalueTreshold <- rbind(PvalueTreshold,
                               list(rownames(PvalueData)[i], PvalueData[, 2][i]))
      }
    }
    colnames(PvalueTreshold) <- c("Genes", "FDR")
    PvalueTreshold <-
      PvalueTreshold[order(PvalueTreshold$FDR),]  # Sort by FDR column
    Genes <- head(PvalueTreshold$Genes, n <- n_genes_to_keep)  # Keep top n rows
    neudata <- data[rownames(data) %in% Genes,]
  }
  if (Pvalue_md == "LRT") {
    count <- 0
    PvalueTreshold <- data.frame(Genes = character(),
                                Fdr = numeric(),
                                stringsAsFactors = FALSE)
    for (i in seq_len(nrow(PvalueData))) {
      if (PvalueData[, 3][i] <= threshold) {
        count <- count + 1
        PvalueTreshold <- rbind(PvalueTreshold,
                               list(rownames(PvalueData)[i], PvalueData[, 2][i]))
      }
    }
    colnames(PvalueTreshold) <- c("Genes", "FDR")
    PvalueTreshold <-PvalueTreshold[order(PvalueTreshold$FDR), ] # Sort by FDR column
    Genes <- head(PvalueTreshold$Genes, n <- n_genes_to_keep)# Keep top n rows
    neudata <- data[rownames(data) %in% Genes, ]
  }
  return(
    list(
      FilterData <- neudata,
      Important_Features <- Genes,
      Statistical_Analysis <- PvalueData)
  )}
