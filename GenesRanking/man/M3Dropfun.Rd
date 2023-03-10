% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Support_VariableGenesMethods.R
\name{M3Dropfun}
\alias{M3Dropfun}
\title{Perform feature selection using M3Drop}
\usage{
M3Dropfun(data, Mt_Method, Mt_threshold)
}
\arguments{
\item{data}{A matrix of count data.}

\item{Mt_Method}{The method to use for selecting mitochondrial genes. Must be one of "M3Drop", "percent", or "counts".}

\item{Mt_threshold}{The threshold to use for selecting mitochondrial genes. If Mt_Method is "M3Drop", this should be a numeric value between 0 and 1. If Mt_Method is "percent", this should be a percentage value between 0 and 100. If Mt_Method is "counts", this should be an integer value.}
}
\value{
A list containing the filtered count data and a list of important features selected by M3Drop.
}
\description{
This function performs feature selection on a matrix of count data using M3Drop.
M3Drop is a method for detecting and filtering out cells with mitochondrial gene expression
patterns, which may be indicative of mitochondrial dysfunction or cell stress.
}
