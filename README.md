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
  - [Use Cases](#use-cases)
  - [Conclusion](#conclusion)



## Introduction
 `GenesRanking` is an R package for filtering and selecting important genes from single-cell RNA-seq data. The package contains two main functions: `Filter_FS_Methods` and `ML_FS_Methods`.

 Our package offers a comprehensive set of tools to facilitate users' first
experience with single-cell RNA sequencing data analysis through feature
selection techniques. The ExampleDataset, Labels, and FilterData objects
are provided to assist users. Labels are character or factor vectors that
identify each cell in the data matrix, while ExampleDataset is a data.frame
object capturing gene expression data for individual cells. FilterData is
a list object containing filtered and/or processed gene expression data.

 To analyze the data, users may choose from two key approaches:
machine learning and statistical methods. However, before applying
any analysis techniques, data preprocessing steps such as normalization 
and/or noise reduction may be required to improve data quality, 
reliability, and to reduce computational time.



## Installation
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(GenesRanking)

```

```{r sessionInfo, warning=FALSE, message=FALSE,results='hide'}
library(GenesRanking)
sessionInfo()
```

## Main Functions

### Filter_FS_Methods


`Filter_FS_Methods` applies Filter feature selection methods. It includes parameters like `HighVariableFilter` and `n_features` for filtering high variability genes, which are often the most informative genes in scRNA-seq data analysis. The `LogTransformation` parameter allows for log transformation of the data using the Seurat package.


Here's an example usage:


```{r Filter_FS_Methods, warning=FALSE, message=FALSE}
# load example data
data(ExampleDataset)
data(Labels)

# apply Filter_FS_Methods
FS_filtered_data <- Filter_FS_Methods(ExampleDataset, Labels, 
                                    n_genes_to_keep=100, 
                                    method="WaldTest",
                                    LogTransformation=TRUE,
                                    HighVariableFIlter=TRUE,
                                    n_features=2000)



```
### ML_FS_Methods

`ML_FS_Methods` applies a specified machine learning method to filter and select important genes from single-cell RNA-seq data.

Here's an example usage:

```{r ML_FS_Methods, warning=FALSE, message=FALSE}
# apply ML_FS_Methods
ML_filtered_data <- ML_FS_Methods(ExampleDataset, Labels, 
                                  n_genes_to_keep=100, 
                                  method="PAM",
                                  LogTransformation=TRUE,
                                  HighVariableFIlter=TRUE,
                                  n_features=1000)

```

## Examination of Results

The results returned by `Filter_FS_Methods` and `ML_FS_Methods` contain a wealth of information that can be leveraged for further analysis. 

```{r Examine Results}

# The names of the list elements provide an overview of the different components of the results:
names(ML_filtered_data)

# For instance, let's take a look at the 'Important_Features' element, which contains the dominant genes identified by the method:
head(ML_filtered_data$Important_Features, 10)

# Another important element could be the 'Feature_Importance_Score', which quantifies the importance of each feature:
head(ML_filtered_data$Feature_Importance_Score, 10)
```


## Use Cases

  The `GenesRanking` package is particularly useful for researchers and bioinformaticians working with single-cell RNA-seq data. By allowing users to focus on the most informative genes and applying specific machine learning methods, it can greatly assist in data analysis and the discovery of new biological insights.
  
## Conclusion

 The `GenesRanking` package provides efficient and user-friendly functions for filtering and selecting important genes from scRNA-seq data. With its capabilities for high variability gene filtering and the application of machine learning methods, it serves as a valuable tool in the analysis of single-cell RNA-seq data.
