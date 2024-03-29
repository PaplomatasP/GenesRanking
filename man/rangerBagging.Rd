% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Support_ML_FIlter.R
\name{rangerBagging}
\alias{rangerBagging}
\title{rangerBagging: Bagging of Ranger models for feature selection}
\usage{
rangerBagging(data, Labels, n_genes_to_keep)
}
\arguments{
\item{data}{A matrix or data frame containing the expression values of genes.}

\item{Labels}{A vector containing the labels corresponding to the samples in data.}

\item{n_genes_to_keep}{An integer representing the number of top important genes to be returned.}
}
\value{
A list containing:
\itemize{
\item FilteredData: A data frame of the original data, filtered to include only the top features.
\item Important_Features: A character vector of the names of the top features.
\item ML_Analysis: A data frame with the importance of each feature.
}
}
\description{
This function applies a bagging ensemble method with the ranger base learner
to a dataset, calculates feature importance and returns the top important features.
}
\examples{
\dontrun{
  results <- rangerBagging(data = data, Labels = labels, n_genes_to_keep = 100)
  head(results$Important_Features)
}
}
