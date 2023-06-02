


#####

methods <- c("LikelihoodRatioTest","WaldTest", "BetaPoissonGLM", "WilcoxonTest",  "KruskalWallis",
             "mRMRe", "QuasiLikelihoodFTest", "KendallCorrelation", "SpearmanCorrelation", "InformationGain",
             "kNN", "svmLinear", "LDA", "RF", "CART",
             "C5.0", "GBM", "xgbTree", "Catboost", "ShapleyValues")
names(filter_data_list_egeod86618)=methods
names(filter_data_list_egeod111727)=methods
names(filter_data_list)=methods
names(filter_data_list_emtab9221)=methods

genes_list1 <- lapply(filter_data_list_egeod86618, function(df) rownames(df))
genes_list2 <- lapply(filter_data_list_egeod111727, function(df) rownames(df))
genes_list3 <- lapply(filter_data_list, function(df) rownames(df)) #results_all_thca_gse_154763
genes_list4 <- lapply(filter_data_list_emtab9221, function(df) rownames(df))
genes_list5 <- lapply(filter_data_list_echad35, function(df) rownames(df))


library(UpSetR)

# Create the upset plot
pl1=upset(fromList(genes_list1), order.by = "freq",nsets = 20,text.scale = 0.7,
          line.size = 0.5,point.size = 1,mb.ratio = c(0.5, 0.5))
pl2=upset(fromList(genes_list2), order.by = "freq",nsets = 20,text.scale = 0.7,
          line.size = 0.5,point.size = 1,mb.ratio = c(0.5, 0.5))
pl3=upset(fromList(genes_list3), order.by = "freq",nsets = 20,text.scale = 0.7,
          line.size = 0.5,point.size = 1,mb.ratio = c(0.5, 0.5))
pl4=upset(fromList(genes_list4), order.by = "freq",nsets = 20,text.scale = 0.7,
          line.size = 0.5,point.size = 1,mb.ratio = c(0.5, 0.5))
pl5=upset(fromList(genes_list5), order.by = "freq",nsets = 20,text.scale = 0.7,
          line.size = 0.5,point.size = 1,mb.ratio = c(0.5, 0.5))


# Convert each plot to grob
pl1_grob <- grid.grabExpr(print(pl1))
pl2_grob <- grid.grabExpr(print(pl2))
pl3_grob <- grid.grabExpr(print(pl3))
pl4_grob <- grid.grabExpr(print(pl4))
pl5_grob <- grid.grabExpr(print(pl5))

# Arrange the plots
grid.arrange(pl1_grob, pl2_grob, pl3_grob, pl4_grob, pl5_grob, ncol = 2, nrow = 3)

