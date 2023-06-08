#' Machine Learning-based Feature Selection
#'
#' This function applies machine learning methods to select top genes based on their importance in classification tasks.
#'
#' @param data A data frame or matrix containing gene expression values. Rows represent genes and columns represent samples.
#' @param Labels A vector containing the labels of the samples. Must be the same length as the number of columns in the data.
#' @param MLmethod Character, the name of the machine learning method to be applied. Supported methods include:
#'   "pam", "gcvEarth", "pls", "rf", "rpart", "C5.0", "gbm", and "xgbTree".
#' @param n_genes_to_keep Integer, the number of top genes to keep after applying the selected machine learning method.
#'
#' @return A list containing the following elements:
#'   - FilterData: A data frame or matrix containing the top n_genes_to_keep genes selected by the specified machine learning method.
#'   - Important_Features: A vector of the names of the top n_genes_to_keep genes.
#'   - ML_Analysis: A data frame containing the importance scores of all genes considered by the machine learning method.
#' @importFrom caret train
#' @importFrom stats predict
ML_filter <- function(data,
                      Labels,
                      MLmethod,
                      n_genes_to_keep) {
    data <- as.data.frame(t(data))
    data$Labels <- as.factor(Labels)

    importanceLimit <- -1
    MLlist <- c(
        "pam",
        "gcvEarth",
        "pls",
        "rf",
        "rpart",
        "C5.0",
        "gbm",
        "xgbTree"
    )
    if (MLmethod %in% MLlist) {
        preProcMethod <- c("knnImpute")
    } else {
        message("ML method input is not valid, exiting...")
        return(-1)
    }
    if (MLmethod == "rpart") {
        colnames(data) <- make.names(colnames(data))
    }

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
        model <-
            train(
                Labels ~ .,
                data = data,
                method = "xgbTree",
                tuneGrid = tune_grid,
                verbose = FALSE,
                metric = "F",
                verbosity = 0
            )
    } else {
        model <- caret::train(
            Labels ~ .,
            data = data,
            method = MLmethod,
            metric = "Accuracy",
            preProc = preProcMethod,
            na.action = stats::na.omit
        )
    }
    message("Model training, please wait!!!!")
    message("Training completed. Searching for the dominant genes!")
    if (MLmethod == "gbm") {
        df_imps <-
            gbm::relative.influence(model$finalModel) # Use gbm::relative.influence for "gbm" method
        df_imps1 <- data.frame(Overall = df_imps)
    } else {
        df_imps <- caret::varImp(model)
        df_imps1 <- df_imps[["importance"]][1]
    }
    df_imps1 <- df_imps1[order(-df_imps1[, 1]), , drop <- FALSE]
    df_imps1 <- subset(df_imps1, df_imps1[, 1] > importanceLimit)
    df_imps2 <- df_imps1[seq_len(n_genes_to_keep), , drop <- FALSE]
    Genes <- rownames(df_imps2)
    Genes <- gsub("`", "", Genes)
    obj <- as.data.frame(t(data))
    newdata <- obj[rownames(obj) %in% Genes, ]
    message("The Analysis is completed.")
    return(list(
        FilterData = newdata,
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
            loss_function = "Logloss",
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

        # Filter the data to keep only the top genes
        filtered_data <- data[, top_genes_indices]
        return(
            list(
                FilteredData = as.data.frame(t(filtered_data)),
                Important_Features = top_genes_indices,
                FeatureImportance = feature_importance
            )
        )
    }


#' ShapleyValues Function
#'
#' This function calculates the Shapley values of a set of features for a given dataset.
#'
#' @param data The input dataset, a matrix or data frame of numeric values.
#' @param Labels A vector of labels for each row in the input data.
#' @param n_genes_to_keep The number of top features to keep.
#'
#' @return A list containing three elements:
#' \describe{
#'   \item{FilterData}{The filtered dataset with only the top features.}
#'   \item{Important_Features}{The names of the top features.}
#'   \item{Statistical_Analysis}{A data frame with the Shapley values for all the features.}
#' }
#'
#' @importFrom ranger ranger
#' @importFrom fastshap explain
ShapleyValuesFun <- function(data, Labels, n_genes_to_keep) {
    rownames(data) <- make.names(rownames(data))

    pfun <-
        function(object, newdata) {
            # extracts prediction from ranger
            predict(object, data = newdata)$predictions[1, ]
        }

    data <- as.data.frame(t(data))

    # Create predictor object
    class_names <- levels(Labels)

    N <- dim(data)[1] # samples

    data$Labels <-
        as.factor(Labels) # add labels to data for model train

    model <-
        ranger::ranger(
            Labels ~ .,
            data = data,
            probability = TRUE,
            num.trees = 100
        )

    # remove labels
    data <- subset(data, select = -Labels)

    # explain contribution of features for each sample
    res <- c()
    for (i in seq_along(class_names)) {
        c <- Labels == class_names[i]
        tmp_res <-
            fastshap::explain(model,
                X = data,
                newdata = data[c, ],
                pred_wrapper = pfun
            )
        res <- rbind(res, tmp_res[i, ]) # keep only correct class
    }
    s_val <- colMeans(res) # mean from all points

    # Initialize rank dataframe
    rank <- data.frame(Genes = colnames(data), Shapley = s_val)

    rank <-
        rank[order(-rank$Shapley), ] # Sort descending (it is already sorted!)

    Genes <- head(rank$Genes, n = n_genes_to_keep) # Keep top n rows

    newdata <- t(data[, colnames(data) %in% Genes])
    newdata <- as.data.frame(newdata)
    return(list(
        FilterData = newdata,
        Important_Features = Genes,
        Statistical_Analysis = rank
    ))
}
