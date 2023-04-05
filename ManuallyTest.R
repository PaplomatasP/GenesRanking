library(GenesRanking)
ls("package:GenesRanking")

#Load availabe data form the GenesRanking Package
data(filtered_data)
data(ExampleDataset)
data(Labels)

#Removes Genes according the chosen Filter
data_filtered <- filter_genes(data=ExampleDataset,filter_method="High_Variable_Genes",n_features= 100)
data_filtered1 <- filter_genes(data=ExampleDataset,filter_method="Remove_Low_Variance",unique_cut= 10)

#Normalization Methods
Norm_data1=normalize_data(ExpressionData=ExampleDataset,normal_method="LogNormalize")
Norm_data2=normalize_data(ExpressionData=ExampleDataset,normal_method="TMM")
Norm_data3=normalize_data(ExpressionData=ExampleDataset,normal_method="CPM")
Norm_data4=normalize_data(ExpressionData=ExampleDataset,normal_method="TPM")
Norm_data5=normalize_data(ExpressionData=ExampleDataset,normal_method="RLE")


#Statistical FS methods
filtered_data1 <- Statistical_Filter(ExampleDataset[1:100,], Pvalue_md = "Waldtest", Labels = Labels,  n_genes_to_keep = 20)
filtered_data2 <- Statistical_Filter(ExampleDataset[1:100,], Pvalue_md = "BPglm", Labels = Labels,  n_genes_to_keep = 20)
filtered_data3 <- Statistical_Filter(ExampleDataset[1:100,], Pvalue_md = "WilcoxonTest", Labels = Labels,  n_genes_to_keep = 20)
filtered_data4 <- Statistical_Filter(ExampleDataset[1:100,], Pvalue_md = "LRT", Labels = Labels,  n_genes_to_keep = 20)

#Non_Tree ML FS methods
Non_Tree1=Non_Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"knn",n_genes_to_keep = 50)
Non_Tree2=Non_Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"lda",n_genes_to_keep = 50)
Non_Tree3=Non_Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"svmRadial",n_genes_to_keep = 50)
Non_Tree4=Non_Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"glmnet",n_genes_to_keep = 50)

#Tree ML FS methods
tree1=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"rf",n_genes_to_keep = 50)
tree2=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"rpart",n_genes_to_keep = 50)
tree3=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"C5.0",n_genes_to_keep = 50)
tree4=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"treebag",n_genes_to_keep = 50)
tree5=Tree_Based_ML(ExampleDataset[1:100,],Labels=Labels,"xgbTree",n_genes_to_keep = 50)


#Ontology Analysis
ontology_results <- Ontology_Analysis(Important_Features = FilterData$Important_Features[1:10],
                                        ontology_type = "GO",
                                         organism = "mouse")



