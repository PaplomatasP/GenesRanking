library(testthat)
library(GenesRanking)

test_that("Statistical_Filter", {
  filtered_data <-
    Statistical_Filter(
      ExampleDataset[1:100, ],
      Pvalue_md = "Waldtest",
      Labels = Labels,
      n_genes_to_keep = 100
    )
  expect_type(filtered_data, "list")
  expect_length(filtered_data, 3)
})

test_that("Non_Tree_Based_ML", {
  result_Wrapper_Based_ML <-
    Non_Tree_Based_ML(
      ExampleDataset[1:50, ],
      Labels = Labels,
      MLmethod = "knn",
      n_genes_to_keep = 15
    )
  expect_type(result_Wrapper_Based_ML, "list")
  expect_length(result_Wrapper_Based_ML, 3)
})

test_that("Tree_Based_ML", {
  result_Tree_Based_ML <-
    Tree_Based_ML(
      ExampleDataset[1:50, ],
      Labels = Labels,
      MLmethod = "C5.0",
      n_genes_to_keep = 15
    )
  expect_type(result_Tree_Based_ML, "list")
  expect_length(result_Tree_Based_ML, 3)
})

test_that("Ontology_Analysis", {
  ontology_results <- Ontology_Analysis(
    Important_Features =
      FilterData$Important_Features[1:10],
    ontology_type = "GO",
    organism = "mouse")
  expect_type(ontology_results, "list")
  expect_length(ontology_results, 10)
})

test_that("filter_genes", {
  data_filtered = filter_genes(ExampleDataset[1:500, 1:100],
                               filter_method = "High_Variable_Genes",
                               n_features = 100)
  expect_type(data_filtered, "list")
  expect_length(data_filtered, 100)
})

test_that("normalize_data", {
  Norm_data = normalize_data(ExpressionData = ExampleDataset,
                             normal_method = "LogNormalize")
  expect_type(Norm_data, "double")
  expect_length(Norm_data, 768000)
})
