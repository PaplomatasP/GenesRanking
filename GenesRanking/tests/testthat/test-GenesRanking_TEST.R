library(testthat)
library(GenesRanking)

#Load availabe data form the GenesRanking Package
data(FilterData)
data(ExampleDataset)
data(Labels)


test_that("Filter_FS_Methods", {
  FS_filtered_data <-
    Filter_FS_Methods(ExampleDataset, Labels, n_genes_to_keep=100,
                     method="WaldTest",
                     LogTransformation=TRUE,
                     HighVariableFIlter=TRUE,
                     n_features=2000)
  expect_type(FS_filtered_data, "list")
  expect_length(FS_filtered_data, 3)
})


test_that("ML_FS_Methods", {
  ML_filtered_data <-
    ML_FS_Methods(ExampleDataset, Labels, n_genes_to_keep=100,
                  method="PAM",
                  LogTransformation=TRUE,
                  HighVariableFIlter=TRUE,
                  n_features=2000)
  expect_type(ML_filtered_data, "list")
  expect_length(ML_filtered_data, 3)
})

