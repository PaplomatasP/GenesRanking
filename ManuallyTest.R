library(BioMarkerX)
ls("package:BioMarkerX")

SCMarker <- Variable_Filter(ExampleDataset[1:500,1:100],  method = "SCMarker")
SelfE <- Variable_Filter(ExampleDataset[1:500,1:100],  method = "SelfE")
DUBStepR <- Variable_Filter(ExampleDataset[1:500,1:100],  method = "DUBStepR")
scPNMF <- Variable_Filter(ExampleDataset[1:500,1:100],  method = "scPNMF")
M3Drop <- Variable_Filter(ExampleDataset[1:500,1:100],  method = "M3Drop")


filtered_data1 <- Statistical_Filter(ExampleDataset[1:100,], Pvalue_md = "Waldtest", Labels = Labels, threshold = 0.05, n_genes_to_keep = 50)
filtered_data2 <- Statistical_Filter(ExampleDataset[1:100,], Pvalue_md = "BPglm", Labels = Labels, threshold = 0.05, n_genes_to_keep = 100)
filtered_data3 <- Statistical_Filter(ExampleDataset[1:100,], Pvalue_md = "WilcoxonTest", Labels = Labels, threshold = 0.05, n_genes_to_keep = 100)
filtered_data4 <- Statistical_Filter(ExampleDataset[1:100,], Pvalue_md = "LRT", Labels = Labels, threshold = 0.05, n_genes_to_keep = 100)

Wrapper1=Wrapper_Based_ML(ExampleDataset[1:100,],Labels=Labels,"knn",importanceLimit = 10,n_genes_to_keep = 50)
Wrapper2=Wrapper_Based_ML(ExampleDataset[1:100,],Labels=Labels,"lda",importanceLimit = 10,n_genes_to_keep = 50)
Wrapper3=Wrapper_Based_ML(ExampleDataset[1:100,],Labels=Labels,"svmRadial",importanceLimit = 10,n_genes_to_keep = 50)
Wrapper4=Wrapper_Based_ML(ExampleDataset[1:100,],Labels=Labels,"glmnet",importanceLimit = 10,n_genes_to_keep = 50)

tree1=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"rf",importanceLimit = 10,n_genes_to_keep = 50)
tree2=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"rpart",importanceLimit = 10,n_genes_to_keep = 50)
tree3=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"C5.0",importanceLimit = 10,n_genes_to_keep = 50)
tree4=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"bartMachine",importanceLimit = 10,n_genes_to_keep = 50)
tree5=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"treebag",importanceLimit = 10,n_genes_to_keep = 50)
tree6=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"xgbTree",importanceLimit = 10,n_genes_to_keep = 50)

Norm_data1=normalize_data(ExpressionData=ExampleDataset,normal_method="LogNormalize")
Norm_data2=normalize_data(ExpressionData=ExampleDataset,normal_method="TMM")
Norm_data3=normalize_data(ExpressionData=ExampleDataset,normal_method="CPM")
Norm_data4=normalize_data(ExpressionData=ExampleDataset,normal_method="TPM")
Norm_data5=normalize_data(ExpressionData=ExampleDataset,normal_method="RLE")


OntologyAnalysis=Ontology_Analysis(FilterData$Important_Features, "KEGG_2021_Human")
library(BioMarkerX)
library(enrichR)

enrichR::enrichr(FilterData$Important_Features,databases ="KEGG_2021_Human")
cdata_filtered <- filter_genes(data=ExampleDataset,filter_method="VariableFeatures",n_features= 100)
data_filtered1 <- filter_genes(data=ExampleDataset,filter_method="nearZeroVar",unique_cut= 10)

