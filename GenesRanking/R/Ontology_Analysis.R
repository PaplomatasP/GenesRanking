#' Perform ontology analysis using clusterProfiler
#'
#' This function performs ontology analysis using clusterProfiler. The user provides a set of important features (e.g., genes), the ontology type to be used (e.g., "KEGG", "GO", "MKEGG"), and the organism (e.g., "human" or "mouse").
#'
#' @param Important_Features A character vector containing the names of the important features (e.g., genes)
#' @param ontology_type A character string indicating the ontology type to be used. It can be one of the following: "KEGG" (default), "GO", "MKEGG".
#' @param organism A character string indicating the organism to be used. It can be either "human" (default) or "mouse".
#' @return A data frame containing the enriched terms and associated statistics
#' @importFrom AnnotationDbi mapIds
#' @importFrom clusterProfiler enrichKEGG enrichGO enrichMKEGG
#' @import org.Hs.eg.db org.Mm.eg.db
#' @export
#' @examples
#' \donttest{
#' # Example with human KEGG pathways
#' data(FilterData)
#' ontology_results <- Ontology_Analysis(Important_Features = FilterData$Important_Features,
#'                                       ontology_type = "KEGG",
#'                                       organism = "human")
#' # Example with mouse Gene Ontology (GO) terms
#' ontology_results <- Ontology_Analysis(Important_Features = FilterData$Important_Features,
#'                                       ontology_type = "GO",
#'                                       organism = "mouse")
#' # Example with mouse KEGG Module (MKEGG) terms
#' ontology_results <- Ontology_Analysis(Important_Features = FilterData$Important_Features,
#'                                       ontology_type = "MKEGG",
#'                                       organism = "mouse")

#' }
Ontology_Analysis <- function(Important_Features, ontology_type = "KEGG", organism = "human") {
  if (organism == "human") {
    org_db <- org.Hs.eg.db
  } else if (organism == "mouse") {
    org_db <- org.Mm.eg.db
  } else {
    stop("Invalid organism specified. Use 'human' or 'mouse'.")
  }

  # Convert gene symbols to Entrez IDs
  entrez_ids <- AnnotationDbi::mapIds(org_db,
                                      keys = Important_Features,
                                      column = "ENTREZID",
                                      keytype = "SYMBOL",
                                      multiVals = "first")

  # Remove NAs
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]

  # Perform enrichment analysis based on the selected ontology type
  if (ontology_type == "KEGG") {
    results <- clusterProfiler::enrichKEGG(gene = entrez_ids,
                                           organism = ifelse(organism == "human", "hsa", "mmu"),
                                           pvalueCutoff = 0.05)
  } else if (ontology_type == "GO") {
    results <- clusterProfiler::enrichGO(gene = entrez_ids,
                                         OrgDb = org_db,
                                         ont = "ALL",
                                         pvalueCutoff = 0.05)
  }  else if (ontology_type == "MKEGG") {
    results <- clusterProfiler::enrichMKEGG(gene = entrez_ids,
                                            organism = ifelse(organism == "human", "hsa", "mmu"),
                                            pvalueCutoff = 0.05)
  } else {
    stop("Invalid ontology_type specified. Use 'KEGG', 'GO', or 'MKEGG'.")
  }
  results=results@result
  return(results)
}
