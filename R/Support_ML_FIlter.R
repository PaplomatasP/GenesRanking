train_model <- function(data, Labels, MLmethod) {
  if (MLmethod == "xgbTree") {
    tune_grid <- expand.grid(
      nrounds = c(50, 100, 150),
      max_depth = c(1, 3, 5),
      eta = c(0.1, 0.05),
      gamma = 0,
      colsample_bytree = 1,
      min_child_weight = 1,
      subsample = 1
    )
  model <- train(
  Labels ~ .,
  data = data,
  method = "xgbTree",
  tuneGrid = tune_grid,
  verbose = FALSE,
  metric = "Accuracy", # Changed from "F" to "Accuracy"
  verbosity = 0
)

  } else if (MLmethod == "rf") {
    model <- caret::train(
      Labels ~ .,
      data = data,
      method = "rf",
      tuneGrid = data.frame(.mtry = sqrt(ncol(data)-1)) # manually specify mtry for rf 
      , ntree = 100
    )
  } else if (MLmethod == "glmnet") {
    
    tune_grid <- expand.grid(alpha = 0, # Ridge regression
                             lambda = seq(0.1, 1, by = 0.1)) # regularization parameter
    
    model <- caret::train(
      Labels ~ .,
      data = data,
      method = "glmnet",
      tuneGrid = tune_grid
    )
  } else if (MLmethod == "C5.0") {
    
    # Set up the parameter grid
    tune_grid <- expand.grid(model = "tree", # tree or rules
                             trials = seq(1, 30, 2), # number of boosting iterations
                             winnow = FALSE ) # whether to use winnowing
   
    model <- caret::train(
      Labels ~ .,
      data = data,
      method = "C5.0",
      tuneGrid = tune_grid
    )
  } else {
    model <- caret::train(
      Labels ~ .,
      data = data,
      method = MLmethod
    )
  }
  return(model)
}



# Function to compute variable importance
compute_varImp <- function(model, MLmethod) {
  if (MLmethod == "gbm") {
    df_imps <- gbm::relative.influence(model$finalModel) # Use gbm::relative.influence for "gbm" method
    df_imps1 <- data.frame(Overall = df_imps)
  } else {
    df_imps <- caret::varImp(model)
    df_imps1 <- df_imps[["importance"]][1]
  }
  
  # Order df_imps1 from larger to smaller
  df_imps1 <- df_imps1[order(-df_imps1[, 1]), , drop = FALSE]
  
  return(df_imps1)
}


# Function to filter important features
filter_important_features <- function(df_imps1, n_genes_to_keep) {
  importanceLimit <- -1
  df_imps1 <- df_imps1[order(-df_imps1[, 1]), , drop <- FALSE]
  df_imps1 <- subset(df_imps1, df_imps1[, 1] > importanceLimit)
  df_imps2 <- df_imps1[seq_len(n_genes_to_keep), , drop <- FALSE]
  Genes <- rownames(df_imps2)
  Genes <- gsub("`", "", Genes)
  return(Genes)
}

#' Machine Learning-based Feature Selection
#'
#' This function applies machine learning methods to select top genes based on their importance in classification tasks.
#'
#' @param data A data frame or matrix containing gene expression values. Rows represent genes and columns represent samples.
#' @param Labels A vector containing the labels of the samples. Must be the same length as the number of columns in the data.
#' @param MLmethod Character, the name of the machine learning method to be applied. Supported methods include:
#'   "pam", "pls", "rf", "rpart", "C5.0", "gbm", "xgbTree" and "glmnet"
#' @param n_genes_to_keep Integer, the number of top genes to keep after applying the selected machine learning method.
#'
#' @return A list containing the following elements:
#'   - FilterData: A data frame or matrix containing the top n_genes_to_keep genes selected by the specified machine learning method.
#'   - Important_Features: A vector of the names of the top n_genes_to_keep genes.
#'   - ML_Analysis: A data frame containing the importance scores of all genes considered by the machine learning method.
#' @importFrom caret train
#' @importFrom stats predict
# Function to train the model
ML_filter <- function(data,
                      Labels,
                      MLmethod,
                      n_genes_to_keep) {
  data <- as.data.frame(t(data))
  data$Labels <- as.factor(Labels)
  
  MLlist <- c(
    "pam",
    "pls",
    "rf",
    "rpart",
    "C5.0",
    "gbm",
    "xgbTree",
    "glmnet"
  )

  # Check if MLmethod is in MLlist
  if (!(MLmethod %in% MLlist)) {
    message("The specified ML method is not supported. Please use one of the following: ", paste(MLlist, collapse = ", "))
    return(-1)
  }
  if (MLmethod == "rpart") {
    colnames(data) <- make.names(colnames(data))
  }

    # Register parallel backend

  #num_cores <- parallel::detectCores() - 4
  #cl <- parallel::makeCluster(num_cores)
  #doParallel::registerDoParallel(cl)
  
  model <- train_model(data, Labels, MLmethod)
  
  #parallel::stopCluster(cl)

  message("Model training, please wait!!!!")
  message("Training completed. Searching for the dominant genes!")
  
  df_imps1 <- compute_varImp(model, MLmethod)
  
  Genes <- filter_important_features(df_imps1, n_genes_to_keep)
  
  obj <- as.data.frame(t(data))
  newdata <- obj[rownames(obj) %in% Genes, ]
  
  newdata <-newdata %>% dplyr::mutate_all(as.numeric)
  
  message("The Analysis is completed.")
  return(list(
    FilteredData  = newdata,
    Important_Features = Genes,
    ML_Analysis = df_imps1
  ))
}


#' Feature selection using CatBoost algorithm
#'
#' This function performs feature selection on single-cell RNA-seq data using the CatBoost algorithm.
#' It trains a CatBoost model, calculates feature importance, and returns the top n_genes_to_keep genes.
#'
#' @param data A matrix or data frame of gene expression data (rows are genes, columns are samples).
#' @param Labels A vector of sample labels (factors or integers) corresponding to the columns in the data.
#' @param n_genes_to_keep The number of top genes to keep based on feature importance.
#'
#' @return A list containing:
#'  \itemize{
#'    \item FilteredData: A data frame or matrix containing only the top n_genes_to_keep genes.
#'    \item Important_Features: A vector of the top n_genes_to_keep gene names.
#'    \item FeatureImportance: A data frame containing the feature importance scores for all genes.
#'  }
#'
#'
CatBoost_feature_selection <-
  function(data, Labels, n_genes_to_keep) {
    data <- as.data.frame(t(data))
    # Convert Labels to integers
    Labels <- as.integer(factor(Labels))
    # Split data into training and testing sets
    train_index <-
      caret::createDataPartition(Labels, p = 0.8, list = FALSE)
    train_data <- data[train_index, ]
    test_data <- data[-train_index, ]
    train_labels <- Labels[train_index]
    test_labels <- Labels[-train_index]
    
    # Create a CatBoost pool
    train_pool <-
      catboost::catboost.load_pool(data = train_data, label = train_labels)
    test_pool <-
      catboost::catboost.load_pool(data = test_data, label = test_labels)
    
    # Train the CatBoost model
    params <- list(
      iterations = 100,
      depth = 6,
      learning_rate = 0.1,
      loss_function = "MultiClass",
      eval_metric = "Accuracy",
      random_seed = 42
    )
    
    model <- catboost::catboost.train(train_pool, test_pool, params)
    # Calculate feature importance
    feature_importance <-
      as.data.frame(
        catboost::catboost.get_feature_importance(model, train_pool,
                                                  type = "FeatureImportance"
        )
      )
    feature_importance <-
      feature_importance[order(-feature_importance[, 1]), , drop <-
                           FALSE]
    # Get the top n_genes_to_keep genes
    top_genes_indices <-
      unlist(head(rownames(feature_importance), n_genes_to_keep))
    colnames(feature_importance)="catboost Score"
    # Filter the data to keep only the top genes
    filtered_data <- data[, top_genes_indices]
    return(
      list(
        FilteredData  = as.data.frame(t(filtered_data)),
        Important_Features = top_genes_indices,
        ML_Analysis = feature_importance
      )
    )
  }



#' rangerBagging: Bagging of Ranger models for feature selection
#'
#' This function applies a bagging ensemble method with the ranger base learner 
#' to a dataset, calculates feature importance and returns the top important features.
#'
#' @param data A matrix or data frame containing the expression values of genes.
#' @param Labels A vector containing the labels corresponding to the samples in data.
#' @param n_genes_to_keep An integer representing the number of top important genes to be returned.
#'
#' @return A list containing:
#'   - FilteredData: A data frame of the original data, filtered to include only the top features.
#'   - Important_Features: A character vector of the names of the top features.
#'   - ML_Analysis: A data frame with the importance of each feature.
#'
#' @importFrom magrittr %>%
#' @importFrom caret train varImp
#' @importFrom ipred bagging
#' @importFrom ranger ranger
#' @examples
#' \dontrun{
#'   results <- rangerBagging(data = data, Labels = labels, n_genes_to_keep = 100)
#'   head(results$Important_Features)
#' }
rangerBagging <- function(data, Labels, n_genes_to_keep) {
  
  data <- as.data.frame(t(data))
  data$Labels <- as.factor(Labels) # add labels to data for model train
  
  # Define the base learner (ranger)
  base_learner <- train(
    x = data ,
    y = Labels,
    method = "ranger"
  )
  
  # Use bagging with ranger as the base learner
  bag <- ipred::bagging(
    formula = Labels ~ .,
    data = data,
    nbagg = 100,
    coob = TRUE,
    fit = base_learner$fit,
    predict = base_learner$predict,
    aggregate = function(data) {
      table(data) %>% prop.table()
    }
  )
  
  importance <- data.frame(imp = caret::varImp(bag))
  genes=rownames(importance)
  # Sort by importance score in descending order
  importance <- data.frame(importance[order(-importance$Overall), ],row.names = genes)
  
  # Select top n_genes_to_keep features
  top_features <- head(rownames(importance), n_genes_to_keep)
  top_features <- gsub("`", '', top_features)


  # Extract the indices of top_features in the data's column names
  indices <- match(top_features, colnames(data))
         
  # Remove NA values if there are any
  indices <- indices[!is.na(indices)]
  
  # Filter the data to keep only the top genes
  filtered_data <- data[, indices] 
  
  return(
    list(
      FilteredData  = as.data.frame(t(filtered_data)),
      Important_Features = top_features,
      ML_Analysis = importance
    )
  )
}
