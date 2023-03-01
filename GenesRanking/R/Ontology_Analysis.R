#' Perform ontology analysis using Enrichr
#'
#' This function performs ontology analysis using Enrichr. The user provides a set of important features (e.g. genes) and the name of the ontology database to be used. The function returns a data frame with the enriched terms and associated statistics.
#'
#' @param Important_Features A character vector containing the names of the important features (e.g. genes)
#' @param Ontology_Terms A character string indicating the name of the ontology database to be used (e.g. "KEGG_2021_Human")
#' @return A data frame containing the enriched terms and associated statistics
#' @importFrom enrichR enrichr
#' @import dplyr
#' @import magrittr
#' @import tidyr
#' @export
#' @examples
#' \donttest{
#' data(FilterData)
#' OntologyAnalysis=Ontology_Analysis(FilterData$Important_Features, "KEGG_2021_Human")
#' }
Ontology_Analysis=function(Important_Features,Ontology_Terms){
   requireNamespace("enrichR", quietly = TRUE)
   options(websiteLive = TRUE)
  dbs <-
    c(
      "KEGG_2021_Human",
      "WikiPathway_2021_Human",
      "BioPlanet_2019",
      "BioCarta_2016",
      "Reactome_2016",
      "MSigDB_Hallmark_2020",
      "GO_Biological_Process_2021",
      "GO_Molecular_Function_2021",
      "GO_Cellular_Component_2021",
      "MGI_Mammalian_Phenotype_Level_4_2021",
      "Human_Phenotype_Ontology",
      "Jensen_DISEASES",
      "DisGeNET",
      "DSigDB",
      "DrugMatrix",
      "OMIM_Disease",
      "HDSigDB_Human_2021",
      "COVID-19_Related_Gene_Sets_2021"
    )


  enriched <- enrichR::enrichr(Important_Features,Ontology_Terms)
  Ontology_Analysis=as.data.frame(enriched)
  return(Ontology_Analysis)
}



