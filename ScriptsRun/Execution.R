options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx32g"))#Restart your R session
library(GenesRanking)
#################### Filter Methods ####################
executionAll_CalcTime_Filter_FS=function(data,Labels,n_genes_to_keep,n_features,
                                       LogTransformation=TRUE,HighVariableFIlter=TRUE){

  # Initialize a data frame to store the method names and execution times

  execution_times_df <- data.frame(Method = rep(NA_character_, 10), Time = rep(NA_real_, 10))


  methods <- c("LikelihoodRatioTest", "BetaPoissonGLM", "WilcoxonTest", "WaldTest", "KruskalWallis",
               "mRMRe", "QuasiLikelihoodFTest", "KendallCorrelation", "SpearmanCorrelation", "InformationGain")
  execution_times_df$Method=methods
 
start_time <- Sys.time()

Filter_FS_Methods1 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="LikelihoodRatioTest"
                                         ,LogTransformation=LogTransformation,
                                         HighVariableFIlter=HighVariableFIlter,n_features=n_features)
 end_time <- Sys.time()

# Calculate the time difference
time_taken <- end_time - start_time
LikelihoodRatioTest_time <- as.numeric(time_taken, units = "secs")
execution_times_df[1,2]=LikelihoodRatioTest_time
start_time <- Sys.time()
Filter_FS_Methods2 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="BetaPoissonGLM"
                                        ,LogTransformation=LogTransformation,
                                        HighVariableFIlter=HighVariableFIlter,n_features=n_features)
end_time <- Sys.time()
# Calculate the time difference
time_taken <- end_time - start_time
BetaPoissonGLM_time <- as.numeric(time_taken, units = "secs")
execution_times_df[2,2]=BetaPoissonGLM_time

start_time <- Sys.time()
Filter_FS_Methods3 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="WilcoxonTest"
                                        ,LogTransformation=LogTransformation,
                                        HighVariableFIlter=HighVariableFIlter,n_features=n_features)
end_time <- Sys.time()
# Calculate the time difference
time_taken <- end_time - start_time
WilcoxonTest_time <- as.numeric(time_taken, units = "secs")
execution_times_df[3,2]=WilcoxonTest_time

start_time <- Sys.time()
Filter_FS_Methods4 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="WaldTest"
                                         ,LogTransformation=LogTransformation,
                                         HighVariableFIlter=HighVariableFIlter,n_features=n_features)
  
end_time <- Sys.time()
# Calculate the time difference
time_taken <- end_time - start_time
WaldTest_time <- as.numeric(time_taken, units = "secs")
execution_times_df[4,2]=WaldTest_time


start_time <- Sys.time()
Filter_FS_Methods5 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="KruskalWallis"
                                        ,LogTransformation=LogTransformation,
                                        HighVariableFIlter=HighVariableFIlter,n_features=n_features)
end_time <- Sys.time()
# Calculate the time difference
time_taken <- end_time - start_time
KruskalWallis_time <- as.numeric(time_taken, units = "secs")
execution_times_df[5,2]=KruskalWallis_time



start_time <- Sys.time()
Filter_FS_Methods6 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="mRMRe"
                                        ,LogTransformation=LogTransformation,
                                        HighVariableFIlter=HighVariableFIlter,n_features=n_features)
end_time <- Sys.time()
# Calculate the time difference
time_taken <- end_time - start_time
mRMRe_time <- as.numeric(time_taken, units = "secs")
execution_times_df[6,2]=mRMRe_time


start_time <- Sys.time()
Filter_FS_Methods7 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="QuasiLikelihoodFTest"
                                        ,LogTransformation=LogTransformation,
                                        HighVariableFIlter=HighVariableFIlter,n_features=n_features)
end_time <- Sys.time()
# Calculate the time difference
time_taken <- end_time - start_time
QuasiLikelihoodFTest_time <- as.numeric(time_taken, units = "secs")
execution_times_df[7,2]=QuasiLikelihoodFTest_time


start_time <- Sys.time()
Filter_FS_Methods8 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="KendallCorrelation"
                                        ,LogTransformation=LogTransformation,
                                        HighVariableFIlter=HighVariableFIlter,n_features=n_features)
end_time <- Sys.time()
# Calculate the time difference
time_taken <- end_time - start_time
KendallCorrelation_time <- as.numeric(time_taken, units = "secs")
execution_times_df[8,2]=KendallCorrelation_time

start_time <- Sys.time()
Filter_FS_Methods9 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="SpearmanCorrelation"
                                        ,LogTransformation=LogTransformation,
                                        HighVariableFIlter=HighVariableFIlter,n_features=n_features)
end_time <- Sys.time()
# Calculate the time difference
time_taken <- end_time - start_time
SpearmanCorrelation_time <- as.numeric(time_taken, units = "secs")
execution_times_df[9,2]=SpearmanCorrelation_time

start_time <- Sys.time()
Filter_FS_Methods10 <<- Filter_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="InformationGain"
                                         ,LogTransformation=LogTransformation,
                                         HighVariableFIlter=HighVariableFIlter,n_features=n_features)
end_time <- Sys.time()
# Calculate the time difference
time_taken <- end_time - start_time
InformationGain_time <- as.numeric(time_taken, units = "secs")
execution_times_df[10,2]=InformationGain_time

return(execution_times_df)
}
# data
# data1=readRDS("ecurd55.rds")
# Labels=readRDS("ecurd55Label.rds")
# Labels=sample(Labels,90000)
# Labels=ifelse(Labels==0,"normal","covid19")
# Labels=as.factor(Labels)
# 
# data1=data1[,1:90000]

Time_FS_Methods=executionAll_CalcTime_Filter_FS(thca_gse_154763,thca_gse_154763_Labels,n_genes_to_keep=100,
                                                n_features=2000,
                                LogTransformation=FALSE,HighVariableFIlter=TRUE)
#saveRDS(Time_FS_Methods,"Time_FS_Methods.rds")


#################### ML Methods ##########################

executionAll_CalcTime_ML_FS=function(data,Labels,n_genes_to_keep,n_features,
                                     LogTransformation=TRUE,HighVariableFIlter=TRUE){

  # Initialize a data frame to store the method names and execution times

  execution_times_df<- data.frame(Method = rep(NA_character_, 10), Time = rep(NA_real_, 10))


  methods <- c("PAM", "gcvEarth","PLS", "RF", "CART",
               "C5.0", "GBM", "xgbTree", "Catboost", "ShapleyValues")
  execution_times_df$Method=methods
  #Statistical FS methods
  start_time <- Sys.time()

  ML_FS_Methods1 <<- ML_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="PAM",
                                   LogTransformation=LogTransformation, HighVariableFIlter=HighVariableFIlter,
                                   n_features=n_features)
  end_time <- Sys.time()

  # Calculate the time difference
  time_taken <- end_time - start_time
  PAM_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[1,2]=PAM_time
  start_time <- Sys.time()
  ML_FS_Methods2 <<- ML_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="gcvEarth",
                                  LogTransformation=LogTransformation,
                                  HighVariableFIlter=HighVariableFIlter, n_features=n_features)
  end_time <- Sys.time()
  # Calculate the time difference
  time_taken <- end_time - start_time
  gcvEarth_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[2,2]=gcvEarth_time

  start_time <- Sys.time()
  ML_FS_Methods3 <<- ML_FS_Methods(data, Labels, n_genes_to_keep=100, method="PLS",
                                  LogTransformation=FALSE,
                                  HighVariableFIlter=TRUE,n_features=2000)
  end_time <- Sys.time()
  # Calculate the time difference
  time_taken <- end_time - start_time
  PLS_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[3,2]=PLS_time

  start_time <- Sys.time()
  ML_FS_Methods4 <<- ML_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="RF",
                                  LogTransformation=LogTransformation,
                                  HighVariableFIlter=HighVariableFIlter,n_features=n_features)
  end_time <- Sys.time()
  # Calculate the time difference
  time_taken <- end_time - start_time
  RF_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[4,2]=RF_time


  start_time <- Sys.time()
  ML_FS_Methods5 <<- ML_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="CART",
                                   LogTransformation=TRUE,
                                  HighVariableFIlter=TRUE,n_features=n_features)
  end_time <- Sys.time()
  # Calculate the time difference
  time_taken <- end_time - start_time
  CART_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[5,2]=CART_time



  start_time <- Sys.time()
  ML_FS_Methods6 <<- ML_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="C5.0",
                                  LogTransformation=LogTransformation,
                                  HighVariableFIlter=HighVariableFIlter,n_features=n_features)
  end_time <- Sys.time()
  # Calculate the time difference
  time_taken <- end_time - start_time
  C5.0_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[6,2]=C5.0_time


  start_time <- Sys.time()
  ML_FS_Methods7 <<- FL
  end_time <- Sys.time()
  # Calculate the time difference
  time_taken <- end_time - start_time
  GBM_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[7,2]=GBM_time


  start_time <- Sys.time()
  ML_FS_Methods8 <<- ML_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="xgbTree",
                                  LogTransformation=LogTransformation,
                                  HighVariableFIlter=HighVariableFIlter,n_features=n_features)
  end_time <- Sys.time()
  # Calculate the time difference
  time_taken <- end_time - start_time
  xgbTree_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[8,2]=xgbTree_time

  start_time <- Sys.time()
  ML_FS_Methods9 <<- ML_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="Catboost",
                                  LogTransformation=LogTransformation,
                                  HighVariableFIlter=HighVariableFIlter,n_features=n_features)
  end_time <- Sys.time()
  # Calculate the time difference
  time_taken <- end_time - start_time
  Catboost_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[9,2]=Catboost_time

  start_time <- Sys.time()
  ML_FS_Methods10 <<- ML_FS_Methods(data, Labels, n_genes_to_keep=n_genes_to_keep, method="ShapleyValues",
                                    LogTransformation=LogTransformation,
                                   HighVariableFIlter=HighVariableFIlter,n_features=n_features)
  end_time <- Sys.time()
  # Calculate the time difference
  time_taken <- end_time - start_time
  ShapleyValues_time <- as.numeric(time_taken, units = "secs")
  execution_times_df[10,2]=ShapleyValues_time

  return(execution_times_df)
}
egeod86618=as.data.frame(t(egeod86618))
Time_ML_Methods=executionAll_CalcTime_ML_FS(egeod86618,egeod86618_Labels,
                                            n_genes_to_keep=100,n_features=2000,
                              LogTransformation=FALSE,HighVariableFIlter=TRUE)


filter_ML_data_list <- c(lapply(list(Filter_FS_Methods1, Filter_FS_Methods2, Filter_FS_Methods3, Filter_FS_Methods4,
                                  Filter_FS_Methods5,Filter_FS_Methods6,Filter_FS_Methods7,
                                  Filter_FS_Methods8,Filter_FS_Methods9,Filter_FS_Methods10,
                                  ML_FS_Methods1,ML_FS_Methods2,ML_FS_Methods3,ML_FS_Methods4,
                                  ML_FS_Methods5,ML_FS_Methods6,ML_FS_Methods7,ML_FS_Methods8,
                                  ML_FS_Methods9,ML_FS_Methods10), function(x) x[[1]]))
Combined_Methods_Time<- rbind(Time_FS_Methods, Time_ML_Methods) 
Combined_Methods_Time$Dataset="thca_gse_154763"
Combined_Methods_Time$Cells=5312

Combined_Methods_Time_echad35=Combined_Methods_Time
filter_data_list_echad35=filter_ML_data_list
names(filter_data_list_echad35)=unlist(Combined_Methods_Time_echad35$Method)


saveRDS(Combined_Methods_Time_echad35,"Combined_Methods_Time_echad35.rds") 
saveRDS(filter_data_list_echad35,"filter_data_list_echad35.rds")
save(list = ls(.GlobalEnv), file = "echad35FS.Rdata")
