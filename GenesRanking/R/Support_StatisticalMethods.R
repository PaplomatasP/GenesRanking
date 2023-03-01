#Function to create SingleCellExperiment  object
#' Function to create a SingleCellExperiment object
#'
#' @param data A matrix or data.frame object representing gene expression data. Rows represent genes, and columns represent samples.
#' @param group A vector of sample labels.
#' @param normal.Method An optional character string indicating the normalization method to use. Default is `NULL`.
#'
#' @return A `SingleCellExperiment` object containing the gene expression data, sample labels, and optional normalized gene expression data.
#'
#' @import SingleCellExperiment
#'
data_fitting = function(data, group, normal.Method) {
  if (missing(normal.Method)) {
    normcounts <- data

  }
  else{
    normcounts <- normalize_data(data, normal.Method)
  }
  SingleCell <-
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = data, normcounts = normcounts))
  gene_df <- data.frame(Gene = rownames(SingleCell))
  cell_df <- data.frame(label = group, cell = colnames(SingleCell))
  rownames(gene_df) <- gene_df$Gene
  rownames(cell_df) <- cell_df$cell
  SingleCell <-
    SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = data, normcounts = normcounts),
      colData = cell_df,
      rowData = gene_df
    )
}

# Function that offers five distinct methods of normalization.


#Calculates the P-value according to the requested method
#' PvalueCalc function
#'
#' This function calculates p-values and FDRs for differential expression analysis using different methods.
#'
#' @param data a SummarizedExperiment object containing gene expression data
#' @param Pvaluemethod a character string specifying the method to be used for calculating p-values and FDRs. Possible values are "Seurat_method", "BPSC_metchod", "DESeq2_method", and "MAST_method".
#' @param Labels a vector of group labels for the samples in the data object
#' @param Test a character string specifying the test to be used in the DESeq2 method. Default is "LRT".
#' @param dmethod a character string specifying the method to be used in the Seurat method. Default is "bimod".
#' @param DESeq2.test a character string specifying the test to be used in the DESeq2 method. Default is "LRT".
#' @param method_MAST a character string specifying the method to be used in the MAST method. Default is "bayesglm".
#' @param MAST.parallel a logical value indicating whether parallel computing should be used in the MAST method. Default is TRUE.
#' @importFrom stats model.matrix na.omit p.adjust
#'
PvalueCalc=function(data,Pvaluemethod,Labels,Test="LRT",dmethod="bimod",DESeq2.test="LRT",method_MAST = "bayesglm",
                    MAST.parallel = TRUE){

  if (Pvaluemethod=="WilcoxonTest"){


    object_Seurat <- as.matrix(SummarizedExperiment::assay(data, "counts"))

    tmp <-Labels

    names(tmp) <- colnames(data)
    meta.data <- data.frame(groups = tmp)

    Seurat.input <- Seurat::CreateSeuratObject(counts = object_Seurat, meta.data = meta.data)

    Seurat.input <- Seurat::NormalizeData(object  = Seurat.input)

    res <- Seurat::FindMarkers(object = Seurat.input, ident.1 = levels(as.factor(tmp))[1],
                               ident.2 = levels(as.factor(tmp))[2], group.by = 'groups',
                               logfc.threshold = -Inf, test.use = dmethod,
                               only.pos = FALSE, verbose = FALSE)
    results_Seurat <- list(gene_names = row.names(res),
                           pvalue = res$p_val,
                           FDR = res$p_val_adj)

    return(results_Seurat)
  }




  else if (Pvaluemethod=="BPglm") {



    object_BPSC <- as.matrix(SummarizedExperiment::assay(data, "normcounts"))

    controllds <- which(Labels == levels(factor(Labels) ) )
    design <- model.matrix(~ Labels)
    resbp <- BPSC::BPglm(data = object_BPSC, controlIds = controllds,
                         design = design, coef = 2, estIntPar = FALSE, useParallel =TRUE )
    FDR <- p.adjust(resbp$PVAL, method = "BH")
    result_BPSC <- list(gene_names = names(resbp$PVAL), pvalue = resbp$PVAL, FDR = FDR)
    return(result_BPSC)
  }



  else if (Pvaluemethod=="LRT") {


    options(warn = -1)
    object_DESeq2 <- as.matrix(SummarizedExperiment::assay(data, "counts"))
    object_DESeq2=na.omit(object_DESeq2)
    object_DESeq2 <- DESeq2::DESeqDataSetFromMatrix(countData = round(object_DESeq2 + 1),
                                                    colData = data.frame(condition = factor(Labels)),
                                                    design = ~ condition)


    if(DESeq2.test == "LRT"){
      object_DESeq2 <- DESeq2::DESeq(object_DESeq2, test = DESeq2.test, parallel = TRUE, betaPrior = FALSE, fitType = "parametric", reduced = ~ 1)
    }
    if (DESeq2.test == "Wald") {
      object_DESeq2 <- DESeq2::DESeq(object_DESeq2, test = DESeq2.test, parallel = TRUE, betaPrior = TRUE, fitType = "parametric")
    }
    # object_DESeq2 <- DESeq2::DESeq(object_DESeq2, test = DESeq2.test, parallel = DESeq2.parallel, ...)
    res_cpm <- DESeq2::results(object_DESeq2)
    result_DESeq2 <- list(gene_names = rownames(res_cpm),
                          pvalue = res_cpm$pvalue,
                          FDR = res_cpm$padj)
    return(result_DESeq2)
  }


  else if (Pvaluemethod=="Waldtest") {



    options(warn = -1)
    object_MAST <- as.matrix(SummarizedExperiment::assay(data, "normcounts"))

    grp <- Labels
    names(grp) <- colnames(object_MAST)

    sca <- MAST::FromMatrix(exprsArray = log2(object_MAST + 1),
                            cData = data.frame(wellKey = names(grp),
                                               grp = grp))
    zlmdata <- MAST::zlm(~ grp, sca, method = method_MAST, parallel = MAST.parallel)
    mast <- MAST::lrTest(zlmdata, "grp")
    FDR  <- p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH")

    result_MAST <- list(gene_names = names(mast[, "hurdle", "Pr(>Chisq)"]),
                        pvalue = mast[, "hurdle", "Pr(>Chisq)"],
                        FDR = FDR)
    return(result_MAST)
  }


}


