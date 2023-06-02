#' Apply a statistical filter to gene expression data for Feature Selection
#' in scRNA-seq Data
#'
#' This function performs a statistical filtering method on gene expression data.
#' It filters the genes based on their statistical significance calculated by
#' various methods.
#'
#' @param data A matrix or data.frame object representing gene expression data.
#'   Rows represent genes, and columns represent Cells.
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
#'
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom BPSC BPglm
#' @importFrom MAST FromMatrix zlm lrTest
#' @importFrom stats model.matrix na.omit p.adjust
#' @importFrom utils head
#' @import Seurat
#'
  Statistical_Filter <- function(data, Pvalue_md, Labels, n_genes_to_keep) {
    options(scipen = 999)
    threshold <- 1
    data1 <- data_fitting(data, Labels)
    PvalueData <- PvalueCalc(data <- data1,
                            Labels <- Labels, Pvaluemethod = Pvalue_md)
    PvalueData  <- as.data.frame(PvalueData)
    PvalueData  <- na.omit(PvalueData)
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
    neudata=as.data.frame(SummarizedExperiment::assay(neudata))
    return(
      list(
        FilterData = neudata,
        Important_Features = Genes,
        Statistical_Analysis = PvalueData)
    )}

  #' Kruskal-Wallis Test for Single-Cell RNA Sequencing Data
  #'
  #' This function performs a Kruskal-Wallis test to rank genes in single-cell RNA sequencing data
  #' based on their expression across different cell types. It returns the filtered data and a summary
  #' of the statistical analysis.
  #'
  #' @param data A numeric matrix or data frame containing the single-cell RNA sequencing data, with
  #'   genes as rows and cells as columns.
  #' @param Labels A factor vector containing the cell type labels for each cell in the data.
  #' @param n_genes_to_keep An integer specifying the number of top-ranked genes to keep in the
  #'   filtered data.
  #'
  #' @return A list containing the following elements:
  #' \itemize{
  #'   \item{FilterData}{A numeric matrix or data frame containing the filtered data with the
  #'     top-ranked genes.}
  #'   \item{Important_Features}{A character vector containing the names of the top-ranked genes.}
  #'   \item{Statistical_Analysis}{A data frame summarizing the Kruskal-Wallis test results, including
  #'     the p-values, false discovery rate (FDR) corrected p-values, and gene names.}
  #' }
  #'
  KruskalWallisFun <- function(data, Labels, n_genes_to_keep) {
    data  <- as.data.frame(t(data))
    data$Labels  <- as.factor(Labels)
    kruskal_w <- apply(data[, -ncol(data)], 2, function(x) kruskal.test(x ~ data$Labels))

    # Extract p-values from the list
    p_values <- sapply(kruskal_w, function(x) x$p.value)

    # FDR Correction
    FDRcor <- (p.adjust(p_values, method = "fdr", n = length(p_values)))

    # Initialize kruskal_w_rank dataframe and add the FDR column
    kruskal_w_rank <- data.frame(P_value = p_values)
    kruskal_w_rank$FDR <- FDRcor
    kruskal_w_rank$Genes <- rownames(kruskal_w_rank)

    colnames(kruskal_w_rank) <- c("P-value", "FDR", "Genes")

    kruskal_w_rank <- kruskal_w_rank[order(kruskal_w_rank$FDR), ] # Sort by FDR column

    Genes <- head(kruskal_w_rank$Genes, n <- n_genes_to_keep) # Keep top n rows

    data <- data[,-ncol(data)]
    data  <- as.data.frame(t(data))
    neudata <- data[rownames(data) %in% Genes, ]

    return(
      list(
        FilterData = neudata,
        Important_Features = Genes,
        Statistical_Analysis = kruskal_w_rank
      )
    )
  }


#' mRMReFun: Minimum Redundancy Maximum Relevance Feature Selection
#'
#' This function performs Minimum Redundancy Maximum Relevance (mRMR) feature selection on a given dataset.
#' It selects the top `n_genes_to_keep` features based on the mRMR criterion.
#'
#' @param data A data frame or matrix containing the expression values with rows as features/genes and columns as samples.
#' @param Labels A character or factor vector containing the labels for each sample.
#' @param n_genes_to_keep An integer specifying the number of top features to select.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{FilterData}{A filtered data frame containing only the selected features/genes.}
#'   \item{Important_Features}{A character vector of the selected feature/genes names.}
#' }
#'
mRMReFun = function(data, Labels, n_genes_to_keep) {
  row.names(data)=make.names(row.names(data))
  data = as.data.frame(t(data))
  data$Labels = as.numeric(as.factor(Labels))

  ## build an mRMRe.Data object
  ge <- mRMRe::mRMR.data(data = data.frame(data))

  # Set the correct target_indices for Labels
  target_index <- ncol(data)

  # Create an mRMR filter and obtain the indices of selected features
  filter <- mRMRe::mRMR.classic("mRMRe.Filter", data = ge, target_indices = target_index,
                                feature_count = n_genes_to_keep)

  ## Get the names of the selected features for each distinct mRMR solution
  Genes <- apply(mRMRe::solutions(filter)[[1]], 2, function(x, y) { return(y[x]) }, y = mRMRe::featureNames(filter))

  ## Filter the data based on the selected genes
  data  <- data[,-ncol(data)]
  data  <- as.data.frame(t(data))
  neudata <- data[rownames(data) %in% Genes, ]
  neudata <- as.data.frame(neudata)
  return(
    list(
      FilterData = neudata,
      Important_Features = Genes

    )
  )
}



#' QLFTfun
#' Filter genes using Quasi-Likelihood F-test (QLFT) from the edgeR package.
#'
#' @param data A SummarizedExperiment object containing gene expression data.
#' @param Labels A vector of character or factor labels for the samples.
#' @param n_genes_to_keep The number of top genes to keep based on the QLFT ranking.
#'
#' @return A list containing:
#' \itemize{
#'   \item{FilterData}{A filtered data matrix with the top n_genes_to_keep genes.}
#'   \item{Important_Features}{A vector of gene names for the top n_genes_to_keep genes.}
#'   \item{Statistical_Analysis}{A data frame containing the genes, P_values, and FDR for the QLFT analysis.}
#' }
QLFTfun <- function(data, Labels, n_genes_to_keep) {
  data <- data_fitting(data, Labels)
  object_edgeR <- as.matrix(SummarizedExperiment::assay(data, "counts"))
  object_edgeR <- data.frame(apply(object_edgeR, 2, function(x) as.numeric(as.character(x))))
  object_edgeR <- na.omit(object_edgeR)

  dge <- edgeR::DGEList(object_edgeR, group = factor(Labels))
  dge <- edgeR::calcNormFactors(dge)
  design <- stats::model.matrix(~Labels)
  dge <- edgeR::estimateDisp(dge, design = design)

  fit <- edgeR::glmQLFit(dge, design = design)
  lrt <- edgeR::glmQLFTest(fit)

  tt <- edgeR::topTags(lrt, n = nrow(object_edgeR))

  result_edgeR <- data.frame(Genes = rownames(tt$table),
                             P_value = tt$table$PValue,
                             FDR = tt$table$FDR)

  result_edgeR <- result_edgeR[order(result_edgeR$FDR), ]

  Genes <- head(result_edgeR$Genes, n = n_genes_to_keep)
  Genes <- as.numeric(Genes)
  Genes  <- rownames(data)[Genes]
  neudata <- data[rownames(data) %in% Genes, ]
  neudata <- as.data.frame(SummarizedExperiment::assay(neudata))
  return(
    list(
      FilterData = neudata,
      Important_Features = Genes,
      Statistical_Analysis = result_edgeR
    )
  )
}
#' Kendall Rank Correlation Feature Selection
#'
#' This function applies Kendall's rank correlation coefficient to select the top n genes
#' with the lowest FDR-adjusted p-values.
#'
#' @param data A data frame or matrix containing gene expression values. Rows represent genes,
#'   and columns represent samples.
#' @param Labels A vector containing the labels of the samples. Must have the same length as
#'   the number of columns in the data.
#' @param n_genes_to_keep The number of top genes to keep based on the FDR-adjusted p-values.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{FilterData}{A data frame or matrix containing the selected genes.}
#'   \item{Important_Features}{A vector of the names of the selected genes.}
#'   \item{Statistical_Analysis}{A data frame containing the p-values, FDR-adjusted p-values, and gene names.}
#' }
#' @importFrom Kendall Kendall
#' @importFrom stats p.adjust
KendallFun <- function(data, Labels, n_genes_to_keep) {


  data = as.data.frame(t(data))
  data$Labels = as.factor(Labels)

  # Calculate Kendall's rank correlation coefficient for each gene
  kendall_result <- apply(data[, -ncol(data)], 2, function(x) Kendall(x, data$Labels))

  # Extract tau and p-values from the list
  kendall_tau <- sapply(kendall_result, function(x) x$tau)
  p_values <- sapply(kendall_result, function(x) x$sl)

  # FDR Correction
  FDRcor <- p.adjust(p_values, method = "fdr", n = length(p_values))

  # Initialize kendall_rank dataframe and add the FDR column
  kendall_rank <- data.frame(P_value = p_values)
  kendall_rank$FDR <- FDRcor
  kendall_rank$Genes <- rownames(kendall_rank)

  colnames(kendall_rank) <- c("P-value", "FDR", "Genes")

  kendall_rank <- kendall_rank[order(kendall_rank$FDR), ] # Sort by FDR column

  Genes <- head(kendall_rank$Genes, n = n_genes_to_keep) # Keep top n rows

  data = data[,-ncol(data)]
  data = as.data.frame(t(data))
  neudata <- data[rownames(data) %in% Genes, ]

  return(
    list(
      FilterData = neudata,
      Important_Features = Genes,
      Statistical_Analysis = kendall_rank
    )
  )
}


#' Correlation-Based Feature Selection
#'
#' This function applies Pearson's correlation coefficient to select the top n genes
#' with the highest absolute correlation coefficients.
#'
#' @param Data A data frame or matrix containing gene expression values. Rows represent genes,
#'   and columns represent samples.
#' @param Labels A vector containing the labels of the samples. Must have the same length as
#'   the number of columns in the Data.
#' @param n_genes_to_keep The number of top genes to keep based on the absolute correlation coefficients.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{FilterData}{A data frame or matrix containing the selected genes.}
#'   \item{Important_Features}{A vector of the names of the selected genes.}
#'   \item{Correlation_Analysis}{A data frame containing the absolute correlation coefficients and gene names.}
#' }
#'
#' @importFrom stats cor
CorrelationFun <- function(Data, Labels, n_genes_to_keep) {
  Data=t(Data)
  Labels=as.numeric(as.factor(Labels))
  # Calculate correlation matrix between features and target variable
  cor_matrix <- cor(Data, Labels, method = "spearman")

  # Extract correlation coefficients with absolute value greater than cutoff
  cor_vals <- abs(cor_matrix[, ncol(cor_matrix)])

  # Sort features by absolute correlation coefficient in descending order
  sorted_features <- sort(cor_vals, decreasing = TRUE)

  # Select top n_genes_to_keep features
  selected_features <- names(sorted_features)[1:n_genes_to_keep]

  # Filter the data based on the selected genes
  filtered_data <- Data[, selected_features]

  return(
    list(
      FilterData = as.data.frame(t(filtered_data)),
      Important_Features = selected_features,
      Correlation_Analysis = as.data.frame(sorted_features)
    )
  )
}




#' Information Gain-Based Feature Selection for Single-Cell RNA-seq Data
#'
#' This function selects the top important features (genes) based on their
#' Information Gain scores from single-cell RNA-seq data.
#'
#' @param Data A data frame or matrix containing the single-cell RNA-seq expression data.
#' @param Labels A vector containing the labels for each sample in the input data.
#' @param n_genes_to_keep An integer specifying the number of top important features (genes) to keep.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{FilterData}{A data frame or matrix containing the filtered input data, keeping only the top important features.}
#'   \item{Important_Features}{A character vector of the top important feature (gene) names.}
#'   \item{Information_Gain_Analysis}{A data frame containing the Information Gain scores for each feature (gene).}
#' }
#'

ig_fun <- function(Data, Labels, n_genes_to_keep) {
  row.names(Data)=make.names(row.names(Data))
  colnames(Data)=make.names(colnames(Data))
  Data<- as.data.frame(t(Data))
  # Combine Data and Labels into a single data frame
  my_dataset <- cbind(Data, Class = Labels)

  # Calculate the Information Gain scores
  ig_scores <- FSelector::information.gain(my_dataset)
  ig_scores$genes=rownames(ig_scores)

  # Sort the features based on their Information Gain scores (in descending order)
  sorted_features <- ig_scores[order(ig_scores$attr_importance, decreasing = TRUE),]

  # Select top n_genes_to_keep features
  selected_features <- rownames(sorted_features)[1:n_genes_to_keep]

  # Filter the data based on the selected features
  filtered_data <- Data[,colnames(Data) %in% selected_features ]
  filtered_data <- as.data.frame(t(filtered_data))

  return(
    list(
      FilterData = filtered_data,
      Important_Features = selected_features,
      Information_Gain_Analysis = sorted_features
    )
  )
}

