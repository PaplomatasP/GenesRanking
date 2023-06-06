library(testthat)
library(GenesRanking)
#Load availabe data form the GenesRanking Package

data(ExampleDataset)
data(Labels)


test_that("ML_filtered_data1", {
  ML_filtered_data <- ML_FS_Methods(ExampleDataset, Labels,
                                    n_genes_to_keep=100,
                                    method="CART",
                                    LogTransformation=FALSE,
                                    HighVariableFIlter=TRUE,
                                    n_features=1000)
  expect_type(ML_filtered_data, "list")
  expect_length(ML_filtered_data, 3)
})


test_that("Filter_FS_Method", {
  FS_filtered_data <- Filter_FS_Methods(ExampleDataset[1:500,], Labels,
                                        n_genes_to_keep=100,
                                        method="mRMRe",
                                        LogTransformation=TRUE,
                                        HighVariableFIlter=TRUE,
                                        n_features=200)
  expect_type(FS_filtered_data, "list")
  expect_length(FS_filtered_data, 2)
})

