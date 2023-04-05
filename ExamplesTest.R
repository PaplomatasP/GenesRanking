library(BioMarkerX)
ls("package:BioMarkerX") 

example("Variable_Filter",package="BioMarkerX")
example("Statistical_Filter",package="BioMarkerX")
example("Wrapper_Based_ML",package="BioMarkerX")
example("Tree_Based_ML",package="BioMarkerX")
example("Ontology_Analysis",package="BioMarkerX")
example("filter_genes",package="BioMarkerX")
example("normalize_data",package="BioMarkerX")


#####################################################
data(ExampleDataset)
filtered_data <- Statistical_Filter(ExampleDataset,  Pvalue_md = "Waldtest", Labels = Labels, threshold = 0.05, n_genes_to_keep = 100)

result <- Tree_Based_ML(ExampleDataset, Labels, MLmethod = "rf", importanceLimit = 10, n_genes_to_keep = 100)

result1 <- Wrapper_Based_ML(ExampleDataset, Labels=Labels, MLmethod = "knn", importanceLimit = 0.05, n_genes_to_keep = 100)


filtered_data <- Statistical_Filter(ExampleDataset, Pvalue_method = "DESeq2_method", Labels = Labels, threshold = 0.01, n_genes_to_keep = 30)

data(filtered_data)
OntologyAnalysis=Ontology_Analysis(filtered_data$Important_Features, "KEGG_2021_Human")
data(ExampleDataset)

result <- Variable_Filter(data=ExampleDataset,  method = "DUBStepR",K=10,Num.pcs=10)
