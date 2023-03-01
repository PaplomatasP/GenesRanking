#' Wrapper_Based_ML for feature selection in scRNA-seq data
#'
#' This function applies machine learning-based feature selection to single-cell RNA sequencing (scRNA-seq) data. The function takes in a data matrix, a vector of cell labels, a machine learning method, and additional parameters to configure the machine learning method. The resulting filtered data and important features are returned as a list.
#'
#' @param data A numeric matrix of gene expression values with rows representing genes and columns representing cells.
#' @param Labels A character or factor vector specifying the label for each cell in \code{data}.
#' @param MLmethod A character string specifying the machine learning method to be used. The valid options are  "knn", "svmRadial",  "glmnet", "lda".
#' @param importanceLimit A numeric value that sets a threshold for the feature importance score. Only features with an importance score greater than this value will be retained.
#' @param n_genes_to_keep An integer value that sets the number of top-ranked features to keep. Only the top \code{n_genes_to_keep} features with the highest importance scores will be retained.
#'
#' @return A list with three elements: \code{FilterData}, \code{Important_Features}, and \code{ML_Analysis}.
#' \code{FilterData} is a numeric matrix of gene expression values with rows representing genes and columns representing cells, after the feature selection process.
#' \code{Important_Features} is a vector of important genes selected by the machine learning-based feature selection method.
#' \code{ML_Analysis} is a data frame containing the importance scores for all features used in the machine learning analysis.
#'
#' @importFrom caret createDataPartition train trainControl varImp
#' @importFrom stats na.omit
#' @examples
#' \donttest{
#' data(ExampleDataset)
#' Wrapper_Based_ML(ExampleDataset,Labels=Labels,"knn",importanceLimit = 10,n_genes_to_keep = 100)
#'}
#' @export
Wrapper_Based_ML = function(data,
                            Labels,
                            MLmethod,
                            importanceLimit,
                            n_genes_to_keep) {

  data=as.data.frame(t(data))
  MLlist = c(

    "knn",
    "svmRadial",
    "glmnet",
    "lda"


  )
  # Set ML parameters, based on the ML method chosen
  if (MLmethod %in% MLlist) {
    preProcMethod <- c("knnImpute")

  } else
  {
    message("ML method input is not valid, exiting...")
    return(-1)
  }


  data$Labels = as.factor(Labels)

  partitionData <-
    caret::createDataPartition(data$Labels, p = 0.75, list = FALSE)
  trainData <- data[partitionData,]
  testData  <- data[-partitionData,]

  trainControl <-
    caret::trainControl(method = "repeatedcv",
                        number = 10,
                        repeats = 3,
                        p = 0.75
    )

  message("Model trainning, please wait!!!!")

  model<-
    caret::train(
      Labels ~ .,
      data = trainData,
      method = MLmethod,  #MLmethod
      metric = "Accuracy",
      preProc = preProcMethod,
      trControl = trainControl,
      na.action = na.omit
    )
  message("Training completed. Searching for the dominant genes!")

  df_imps <- varImp(model)



  df_imps1 <- df_imps[["importance"]][1]
  df_imps1 <- df_imps1[order(-df_imps1[, 1]), , drop = FALSE]
  df_imps2 <- df_imps1[1:n_genes_to_keep, , drop = FALSE]
  Genes=rownames(df_imps2)
  Genes <- gsub("`", "", Genes)

  obj<-as.data.frame(t(data))
  newdata = obj[  rownames(obj) %in% Genes,]

  message("Training completed. Searching for the dominant genes!")

  return(list(FilterData=newdata,Important_Features=Genes,ML_Analysis=df_imps1)  )
}

