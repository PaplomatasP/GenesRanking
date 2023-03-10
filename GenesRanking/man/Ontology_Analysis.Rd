% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ontology_Analysis.R
\name{Ontology_Analysis}
\alias{Ontology_Analysis}
\title{Perform ontology analysis using Enrichr}
\usage{
Ontology_Analysis(Important_Features, Ontology_Terms)
}
\arguments{
\item{Important_Features}{A character vector containing the names of the important features (e.g. genes)}

\item{Ontology_Terms}{A character string indicating the name of the ontology database to be used (e.g. "KEGG_2021_Human")}
}
\value{
A data frame containing the enriched terms and associated statistics
}
\description{
This function performs ontology analysis using Enrichr. The user provides a set of important features (e.g. genes) and the name of the ontology database to be used. The function returns a data frame with the enriched terms and associated statistics.
}
\examples{
\donttest{
data(FilterData)
OntologyAnalysis=Ontology_Analysis(FilterData$Important_Features, "KEGG_2021_Human")
}
}
