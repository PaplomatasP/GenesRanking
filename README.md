---
title: "GenesRanking_vignette"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{GenesRanking_vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```
# GenesRanking

  - [GenesRanking](#GenesRanking)
  - [Introduction](#introduction)
  - [Installation](#installation)
  - [Main Functions](#main-functions)
  - [Examination of Results](#examination-of-results)
  - [Use Cases and Conclusion](#use-cases-and-conclusion)



## Introduction
 
The `GenesRanking` is an R package designed to filter and select significant genes from single-cell RNA-seq data. This package consists of two key functions: `Filter_FS_Methods` and `ML_FS_Methods`.

The GenesRanking package offers a comprehensive set of tools to facilitate the initial experience of users with single-cell RNA sequencing data analysis using feature selection techniques. We provide `ExampleDataset` and `Labels` objects to assist users in the process. Labels are character or factor vectors that identify each cell in the data matrix, while the ExampleDataset is a data.frame object that encapsulates gene expression data for individual cells.

Data analysis in various domains offers two primary approaches: machine learning and statistical methods. Nonetheless, before applying any analysis techniques, it is essential to perform preprocessing steps like normalization and noise reduction. These steps improve the data's quality and reliability while reducing computational time.

In light of this, our functions offer the options for `LogTransformation` and `HighVariableFilter`. These choices enable users to apply log transformation to the data and utilize a high variability filter, respectively.

## Installation
```
install.packages("GenesRanking")

```

```{r sessionInfo, warning=FALSE, message=FALSE,results='hide'}
library(GenesRanking)
library(magrittr)
#Please installe the BPSC and catboost packages as refers below:
#install.packages("devtools")
library("devtools")
#install_github("nghiavtr/BPSC")
#install.packages('devtools')
#devtools::install_url('https://github.com/catboost/catboost/releases/download/v1.2/catboost-R-Windows-1.2.tgz', INSTALL_opts = c("--no-multiarch", "--no-test-load"))

library("catboost")
library("BPSC")

```

## Main Functions

### Filter_FS_Methods

The `Filter_FS_Methods` function applies filter feature selection methods based on a statistical approach. It includes parameters such as HighVariableFilter and n_features to filter high variability genes, which are often the most informative genes in scRNA-seq data analysis. The LogTransformation parameter allows for log transformation of the data. Both of these parameters utilize the Seurat package for their implementation.

Available methods for the Filter_FS_Methods function include `WaldTest, BetaPoissonGLM, WilcoxonTest, LikelihoodRatioTest, KruskalWallis, mRMRe, QuasiLikelihoodFTest, KendallCorrelation, SpearmanCorrelation, and InformationGain`.


Here's an example usage:


```{r Filter_FS_Methods, warning=FALSE, message=FALSE}
# load example data
data(ExampleDataset)
data(Labels)

# apply Filter_FS_Methods
FS_filtered_data <- Filter_FS_Methods(ExampleDataset[1:500,], Labels, 
                                    n_genes_to_keep=100, 
                                    method="InformationGain",
                                    LogTransformation=TRUE,
                                    HighVariableFIlter=TRUE,
                                    n_features=200)



```
### ML_FS_Methods

The `ML_FS_Methods` function applies specified machine learning algorithms to filter and select important genes from single-cell RNA-seq data. This function uses various machine learning techniques to identify genes that are most relevant for a specific task or analysis. By leveraging the power of machine learning algorithms, it aims to provide a data-driven approach for gene selection and feature extraction in scRNA-seq data analysis. Additionally, the function also offers the availability of HighVariableFilter and LogTransformation options.

Available methods for the ML_FS_Methods function include `PAM, glmnet, PLS, RF, CART, C5.0, GBM, xgbTree, Catboost, and rangerBagging`.

Here's an example usage:

```{r ML_FS_Methods, warning=FALSE, message=FALSE}
# apply ML_FS_Methods
ML_filtered_data <- ML_FS_Methods(ExampleDataset[1:500,], Labels, 
                                  n_genes_to_keep=100, 
                                  method="rangerBagging",
                                  LogTransformation=TRUE,
                                  HighVariableFIlter=TRUE,
                                  n_features=200)

```

## Examination of Results

The results returned by `Filter_FS_Methods` and `ML_FS_Methods` contain a wealth of information that can be leveraged for further analysis. 

```{r Examine Results}

# The names of the list elements provide an overview of the different components of the results:
names(ML_filtered_data)

# For instance, let's take a look at the 'Important_Features' element, which contains the dominant genes identified by the ML FS method:
head(ML_filtered_data$Important_Features, 10)



#let's take a look at the 'Important_Features' element, which contains the dominant genes identified by the Filter FS method:
head(FS_filtered_data$Important_Features, 10)

# Another important element could be the 'Feature_Importance_Score', which quantifies the importance of each feature:
head(FS_filtered_data[[3]], 10)
```




## Use Cases and Conclusion

The `GenesRanking package` serves as a practical tool for researchers and bioinformaticians working with gene expression data. This package simplifies the process of identifying the most informative genes and allows users to apply specific machine learning or filter feature selection methods, enhancing both data analysis and the discovery of new biological insights.

GenesRanking provides user-friendly and efficient functions for filtering and selecting important genes from the complex landscape of gene expression data. It offers unique capabilities, such as high variability gene filtering and the application of various machine learning methods. Furthermore, it provides a range of 20 different filter and machine learning techniques, enabling users to effectively identify potential biomarkers and deepen their understanding of biological processes.

In conclusion, GenesRanking is a versatile and powerful tool for anyone looking to navigate the intricate field of gene expression data analysis.
