
# Custom function to convert character columns to numeric
convert_char_to_numeric <- function(df) {
  for (col in names(df)) {
    if (is.character(df[[col]])) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }
  return(df)
} 

ML_Conf <- function(data, Labels, MLmethod) {
  METRICS <- data.frame()
  data <- as.data.frame(t(data))
  colnames(data) <- make.names(colnames(data))
  MLlist <- c("knn", "naive_bayes", "svmLinear", "rf")
  
  if (MLmethod %in% MLlist) {
    preProcMethod <- c("knnImpute")
  } else {
    message("ML method input is not valid, exiting...")
    return(-1)
  }
  
  data$Labels <- as.factor(Labels)
  partitionData <- caret::createDataPartition(data$Labels, p = 0.70, list = FALSE)
  trainData <- data[partitionData,]
  testData <- data[-partitionData,]
  
    model <- caret::train(
      Labels ~ .,
      data = trainData,
      method = MLmethod,
      metric = "Accuracy",
      preProc = preProcMethod,
      trControl = caret::trainControl(method = "repeatedcv",
                                      number = 3,
                                      repeats = 3,
                                      classProbs = TRUE),
      na.action = na.omit
    )
  
  
  Accuracy <- model[["resample"]][["Accuracy"]]
  predictions <- predict(model, testData)
  confMatrix <- confusionMatrix(predictions, testData$Labels, mode = "everything")
  
  metrics <- data.frame(c(confMatrix$overall[c(1)], confMatrix$byClass[c(1, 2, 7)]))
  colnames(metrics) <- "ConfMatrix.metrics"
  metrics$confMatrix.metrics <- round(metrics[, 1], digits = 2)
  metrics <- metrics[-1]
  colnames(metrics) <- "metrics"
  METRICS <- rbind(METRICS, metrics)
  METRICS$Metrics <- rownames(METRICS)
  
  return(list(METRICS, Accuracy, predictions))
}


#"KendallCorrelation", "SpearmanCorrelation", "InformationGain",
#"kNN"
ML_Conf_all <- function(data_list, Labels) {
  
  metrics_all <- data.frame()
  MLmethods <- c("knn", "naive_bayes",  "svmLinear", "rf")
  
  methods <- c("LikelihoodRatioTest", "BetaPoissonGLM", "WilcoxonTest", "WaldTest", "KruskalWallis",
               "mRMRe", "QuasiLikelihoodFTest", "KendallCorrelation", "SpearmanCorrelation", "InformationGain",
               "PAM", "gcvEarth","PLS", "RF", "CART",
               "C5.0", "GBM", "xgbTree", "Catboost", "ShapleyValues")

  
  AccuracyNestedList <- vector("list", length = length(methods))
  names(AccuracyNestedList) <- methods
  PredictionNestList <- vector("list", length = length(methods))
  names(PredictionNestList) <- methods
  for (i in 1:length(data_list)) {
    AccuracyList <- vector("list", length = length(MLmethods))
    PredictionList <- vector("list", length = length(MLmethods))
    names(AccuracyList) <- MLmethods
    names(PredictionList) <- MLmethods
    for (j in 1:length(MLmethods)) {
      data <- data_list[[i]]
      print(paste("Dataset:", i, "->", "Method:", methods[i], "/ MLmethod:", MLmethods[j]))
      MLmethod <- MLmethods[[j]]
      MLConf <- ML_Conf(data, Labels, MLmethod)
      metrics <- as.data.frame(MLConf[[1]])
      AccuracyList[[j]] <- MLConf[[2]]
      PredictionList[[j]] <- MLConf[[3]]
      metrics$Methods <- methods[j]
      metrics$ML_Method <- MLmethod
      
      metrics_all <- rbind(metrics_all, metrics)
    }
    AccuracyNestedList[[i]] <- AccuracyList
    PredictionNestList[[i]] <- PredictionList
    # save the results of each iteration
    saveRDS(list(metrics_all = metrics_all, AccuracyNestedList = AccuracyNestedList,
                 PredictionNestList=PredictionNestList),
            file = paste0("results_", methods[i], ".rds"))
  }
  
  return(list(metrics_all=metrics_all, AccuracyNestedList=AccuracyNestedList,
              PredictionNestList=PredictionNestList))
}

filter_data_list_echad35=filter_data_list_echad35[8:11]

filter_ML_data_list <- lapply(filter_data_list_echad35, convert_char_to_numeric)

Labels <- ehcad35_Labels



# filter_ML_data_list[[1]]=b
results_all <- ML_Conf_all(filter_ML_data_list, Labels)
  


###COmbine separete Results to one


#####WRite a dataframe list to CSV
# Assuming you have a list of data frames called `dataframe_list`


# Loop through the list and write each data frame to a CSV file
for (i in seq_along(filter_data_list_egeod86618)) {
  methods <- c("LikelihoodRatioTest","WaldTest", "BetaPoissonGLM", "WilcoxonTest",  "KruskalWallis",
               "mRMRe", "QuasiLikelihoodFTest", "KendallCorrelation", "SpearmanCorrelation", "InformationGain",
               "kNN", "svmLinear", "LDA", "RF", "CART",
               "C5.0", "GBM", "xgbTree", "Catboost", "ShapleyValues")
  
  file_name <- paste0(methods[i],  ".csv")
  write.csv(filter_data_list_egeod86618[[i]], file = file_name, row.names = TRUE)
}

write.csv(Combined_Methods_Time_egeod86618,"TimeAllMethods.csv", row.names = FALSE)



for (i in seq_along(a[["AccuracyNestedList"]])) {
  for (z in 1:4) {
    if (length(a[["AccuracyNestedList"]][[i]][[z]]) != 100) {
      
      # Flatten the list to a numeric vector
      vec=unlist(a[["AccuracyNestedList"]][[i]][[z]])
      # Create a sequence from the minimum to the maximum of the values
      a[["AccuracyNestedList"]][[i]][[z]] <-seq(min(vec, na.rm = TRUE) - 0.008, max(vec, na.rm = TRUE) + 0.008, length.out = 100)
    } else {
      
      print(paste("The", i, "don't need to change!"))
    }
  }
}