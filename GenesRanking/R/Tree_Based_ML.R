#' Applies tree-based machine learning for feature selection in scRNA-seq data
#'
#' This function applies machine learning-based feature selection using
#' tree-based models to single-cell RNA sequencing (scRNA-seq) data. The function
#' takes as input a data matrix, a vector of cell labels, a machine learning
#' method, and additional parameters to configure the machine learning method.
#' The resulting filtered data, important features, and the ML analysis are
#' returned as a list.
#'
#' @param data A numeric matrix of gene expression values with rows representing
#'   genes and columns representing cells.
#' @param Labels A character or factor vector specifying the label for each cell
#'   in \code{data}.
#' @param MLmethod A character string specifying the machine learning method to
#'   be used. The valid options are "rf", "knn","rpart","C5.0","bartMachine",
#'   "treebag", or "xgbTree".
#' @param n_genes_to_keep An integer value that sets the number of top-ranked
#'   features to keep. Only the top \code{n_genes_to_keep} features with the
#'   highest importance scores will be retained.
#'
#' @return A list with three elements: \code{FilterData},
#'   \code{Important_Features}, and \code{ML_Analysis}.
#' \code{FilterData} is a numeric matrix of gene expression values with rows
#'   representing genes and columns representing cells, after the feature
#'   selection process.
#' \code{Important_Features} is a vector of important genes selected by the
#'   machine learning-based feature selection method.
#' \code{ML_Analysis} is a data frame containing the importance scores for all
#'   features used in the machine learning analysis.
#'
#' @importFrom caret createDataPartition train trainControl varImp
#' @importFrom stats na.omit
#' @examples
#'
#' data(ExampleDataset)
#' Tree_Based_ML(head(ExampleDataset, 10), Labels, MLmethod<-"rf",
#'   n_genes_to_keep<-8)
#'
#' @export
Tree_Based_ML <- function(data,
                          Labels,
                          MLmethod,
                          n_genes_to_keep) {
  importanceLimit <- -1
  data <- as.data.frame(t(data))
  MLlist <- c(
    "rf","rpart","C5.0","treebag","xgbTree"
  )
  if (MLmethod %in% MLlist) {
    preProcMethod <- c("knnImpute")
  } else {
    message("ML method input is not valid, exiting...")
    return(-1)
  }
  if (MLmethod == "rpart" | MLmethod == "treebag") {
    colnames(data) <- make.names(colnames(data))
  }
  data$Labels <- as.factor(Labels)
  partitionData <- caret::createDataPartition(
    data$Labels, p = 0.75, list = FALSE)
  trainData <- data[partitionData,]
  testData  <- data[-partitionData,]
  trainControl <- caret::trainControl(method = "repeatedcv",
                                      number = 2,repeats = 1,
                                      p = 0.75)
  message("Model training, please wait!!!!")
  model <- caret::train(
    Labels ~ .,data <- trainData,
    method = MLmethod,
    metric = "Accuracy",
    preProc = preProcMethod,
    trControl = trainControl,
    na.action = stats::na.omit)
  message("Training completed. Searching for the dominant genes!")
  df_imps <- varImp(model)
  df_imps1 <- df_imps[["importance"]][1]
  df_imps1 <- df_imps1[order(-df_imps1[, 1]), , drop <- FALSE]
  df_imps1 <- subset(df_imps1, df_imps1[, 1] > importanceLimit)
  df_imps2 <- df_imps1[seq_len(n_genes_to_keep), , drop <- FALSE]
  Genes <- rownames(df_imps2)
  Genes <- gsub("`", "", Genes)
  obj <- as.data.frame(t(data))
  newdata <- obj[rownames(obj) %in% Genes,]
  message("The Analysis is completed.")
  return(
    list(
      FilterData <- newdata,
      Important_Features <- Genes,
      ML_Analysis <- df_imps1
    ))
}
