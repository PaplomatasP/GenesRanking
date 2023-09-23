#' ML-Based Feature Selection Methods
#'
#' This function applies machine learning (ML) methods for feature selection in single-cell
#' RNA-seq data. It takes as input a gene expression matrix, sample labels,
#' the number of top genes to retain, and the name of the ML method.
#' Optional parameters include log transformation and high variability gene filtering.
#' The function returns a list containing the filtered data, key features, and other
#' information relevant to the selected method.
#' Available ML methods include PAM, rangerBagging, PLS, RF, CART, C5.0, GBM, xgbTree,
#' Catboost, and glmnet. These methods can help identify the most informative genes
#' for downstream analysis.
#' This function is designed to facilitate the comparison of different ML methods
#' and to assist in the identification of the most informative genes for further analysis.
#'
#' @param data Matrix or data frame of input data (features and samples).
#' @param Labels Vector of true class labels of the samples.
#' @param n_genes_to_keep Integer. Number of top-ranked genes to retain.
#' @param method Character. ML method for feature selection. Supported methods: 
#' "PAM", "rangerBagging", "PLS", "RF", "CART", "C5.0", "GBM", "xgbTree", "Catboost", 
#' "glmnet".
#' @param LogTransformation Logical. If TRUE, apply log transformation using Seurat.
#' Default is TRUE.
#' @param HighVariableFIlter Logical. If TRUE, apply high variability gene filtering. 
#' Default is TRUE.
#' @param n_features Integer. Number of genes with high variability to retain.
#'
#' @return List containing the filtered data, key features, and other information 
#' relevant to the selected method.
#' @examples
#' # Example usage:
#' # Apply the rangerBagging method and keep the top 100 genes
#' result <- ML_FS_Methods(ExampleDataset[1:200, ], Labels,
#'     n_genes_to_keep = 50, method = "rangerBagging",
#'     LogTransformation = TRUE, HighVariableFIlter = TRUE, n_features = 100)
#' @export
ML_FS_Methods <-
    function(data,
             Labels,
             n_genes_to_keep,
             method,
             LogTransformation = TRUE,
             HighVariableFIlter = TRUE,
             n_features) {
  min_genes <- nrow(data)/4
  min_cell_percentage <- 0.1
  keep_going <- TRUE
  
  while (keep_going) {
    # Filter lowly-expressed genes
    min_cells <- ceiling(min_cell_percentage * ncol(data)) # At least min_cell_percentage of cells
    gene_counts <- rowSums(data > 0)
    filtered_data <- data[gene_counts >= min_cells, ]
    
    # Check if we have at least nrow(data)/4 genes
    if (nrow(filtered_data) >= min_genes) {
      keep_going <- FALSE
    } else {
      # If not, reduce min_cell_percentage to keep more genes
      min_cell_percentage <- min_cell_percentage - 0.01
    }
  }
#if (n_features < min_genes){
#      message("The number of genes ")
   # }
        
  # Set data to the filtered data
  data <- filtered_data
        if (HighVariableFIlter == TRUE) {
            seurat_obj <- Seurat::CreateSeuratObject(data)
            data_seurat <- Seurat::FindVariableFeatures(seurat_obj,
                selection.method = "vst",
                nfeatures = n_features
            )
            genes <- head(Seurat::VariableFeatures(data_seurat), n_features)
            data <- as.data.frame(data[rownames(data) %in% genes, ])
        }
        if (LogTransformation == TRUE) {
            scale_factor <- mean(colSums(data))
            data <- as.matrix(Seurat::LogNormalize(data, scale.factor <-
                scale_factor))
        }
        switch(method,
            "PAM" = {
                result <-
                    ML_filter(
                        data,
                        MLmethod = "pam",
                        Labels = Labels,
                        n_genes_to_keep = n_genes_to_keep
                    )
            },
            "rangerBagging" = {
                result <-
                  rangerBagging(
                        data,
                        Labels = Labels,
                        n_genes_to_keep = n_genes_to_keep
                    )
            },
            "PLS" = {
                result <-
                    ML_filter(
                        data,
                        MLmethod = "pls",
                        Labels = Labels,
                        n_genes_to_keep = n_genes_to_keep
                    )
            },
            "RF" = {
                result <-
                    ML_filter(
                        data,
                        MLmethod = "rf",
                        Labels = Labels,
                        n_genes_to_keep = n_genes_to_keep
                    )
            },
            "CART" = {
                result <-
                    ML_filter(
                        data,
                        MLmethod = "rpart",
                        Labels = Labels,
                        n_genes_to_keep = n_genes_to_keep
                    )
            },
            "C5.0" = {
                result <-
                    ML_filter(
                        data,
                        MLmethod = "C5.0",
                        Labels = Labels,
                        n_genes_to_keep = n_genes_to_keep
                    )
            },
            "GBM" = {
                result <-
                    ML_filter(
                        data,
                        MLmethod = "gbm",
                        Labels = Labels,
                        n_genes_to_keep = n_genes_to_keep
                    )
            },
            "xgbTree" = {
                result <-
                    ML_filter(
                        data,
                        MLmethod = "xgbTree",
                        Labels = Labels,
                        n_genes_to_keep = n_genes_to_keep
                    )
            },
            "Catboost" = {
                result <-
                    CatBoost_feature_selection(data, Labels = Labels, n_genes_to_keep = n_genes_to_keep)
            },
            "glmnet" = {
              result <-
                ML_filter(
                  data,
                  MLmethod = "glmnet",
                  Labels = Labels,
                  n_genes_to_keep = n_genes_to_keep
                )
            },
            {
                stop("Invalid method. Please choose a valid feature selection method.")
            }
        )

        return(result)
    }
