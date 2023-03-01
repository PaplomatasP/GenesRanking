library(testthat)
library(GenesRanking)
library(enrichR)

test_that("Variable_Filter", {
  data(ExampleDataset)
  result <- Variable_Filter(ExampleDataset[1:500,1:100],  method = "DUBStepR", K=10,Num.pcs=10)
  expect_type(result,"list")
  expect_length(result,3)


})


test_that("Statistical_Filter", {
  data(ExampleDataset)
  filtered_data <- Statistical_Filter(ExampleDataset[1:100,], Pvalue_md = "Waldtest", Labels = Labels, threshold = 0.05, n_genes_to_keep = 100)
  expect_type(filtered_data,"list")
  expect_length(filtered_data,3)

})

test_that("Wrapper_Based_ML", {
  data(ExampleDataset)
  data(Labels)
  result_Wrapper_Based_ML <- Wrapper_Based_ML(ExampleDataset[1:50,], Labels=Labels, MLmethod = "knn", importanceLimit = 10, n_genes_to_keep = 15)
  expect_type(result_Wrapper_Based_ML,"list")
  expect_length(result_Wrapper_Based_ML,3)
})

test_that("Tree_Based_ML", {
  data(ExampleDataset)
  data(Labels)
  result_Tree_Based_ML <- Tree_Based_ML(ExampleDataset[1:50,], Labels=Labels, MLmethod = "C5.0", importanceLimit = 10, n_genes_to_keep = 15)
  expect_type(result_Tree_Based_ML,"list")
  expect_length(result_Tree_Based_ML,3)
})

test_that("Ontology_Analysis", {
  data(FilterData)
  OntologyAnalysis=Ontology_Analysis(FilterData$Important_Features, "KEGG_2021_Human")
  expect_type(OntologyAnalysis,"list")
  expect_length(OntologyAnalysis,9)
})

test_that("filter_genes", {
  data(ExampleDataset)
  data_filtered <- filter_genes(data=ExampleDataset[1:500,1:100], filter_method = "VariableFeatures", n_features = 100)
  expect_type(data_filtered,"list")
  expect_length(data_filtered,100)
})

test_that("normalize_data", {
  data(ExampleDataset)
  Norm_data=normalize_data(ExpressionData=ExampleDataset, normal_method="LogNormalize")
  expect_type(Norm_data,"double")
  expect_length(Norm_data,768000)
})
