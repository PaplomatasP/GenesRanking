CorMetadata=read.csv("human_sfg_corr_meta.csv",row.names = 1)
Cordata=read.csv("human_sfg_corr.csv",row.names = 1)
Cordata=as.data.frame(t(Cordata))


TableFunctionFS=function(data,Labels,n,Index=FALSE){
  scmap <- Variable_Filter(data,  method = "scmap",labels=Labels,nGenes=n)

  filtered_data <- Statistical_Filter(data, Pvalue_md = "Waldtest", Labels = Labels,  n_genes_to_keep = n)
  #Wrapper=Non_Tree_Based_ML(data,Labels<-Labels,"knn",n_genes_to_keep = n)
  #tree<-Tree_Based_ML(data,Labels=Labels,"C5.0",n_genes_to_keep = n)

  RankDataFrame=data.frame(matrix(nrow=n,ncol=2))
  l1=length(scmap$Important_Features):n
  RankDataFrame$X1[1:l1]=scmap$Important_Features

  l2=length(filtered_data$Important_Features):n
  RankDataFrame$X2[1:l2]=filtered_data$Important_Features

  #l3=length(Wrapper$Important_Features):n
  #RankDataFrame$X3[1:l3]=Wrapper$Important_Features

  # l4=length(tree$Important_Features):n
  #RankDataFrame$X4[1:l4]=tree$Important_Features
  colnames(RankDataFrame)=c("scmap","Waldtest")

  if (Index==TRUE){
    RankDataFrame$scmap= match(unlist(RankDataFrame[,1]),rownames(data))
    RankDataFrame$Waldtest= match(unlist(RankDataFrame[,2]),rownames(data))
    #RankDataFrame[,1]= match(unlist(RankDataFrame[,1]),rownames(data))
    #RankDataFrame[,1]= match(unlist(RankDataFrame[,1]),rownames(data))

  }
  return(RankDataFrame)


}


RankCorData=TableFunctionFS(Cordata,Labels=as.factor(CorMetadata$group) ,n=2000,Index=TRUE)


RankCorData=TableFunctionFS(ExampleDataset[1:500,],Labels=Labels ,n=500,Index=TRUE)

scmap <- Variable_Filter(Cordata,  method = "scmap",labels=as.factor(CorMetadata$group),nGenes=2000)

Cordata=TableFunctionFS(data=Cordata,Labels=CorMetadata$group,n=2000)


scPNMF <- Variable_Filter(ExampleDataset[1:500,],DM="EucDist",nGenes=500,  method = "scPNMF")
filtered_data <- Statistical_Filter(ExampleDataset[1:500,], Pvalue_md = "Waldtest", Labels = Labels, n_genes_to_keep = 500)
Wrapper=Non_Tree_Based_ML(ExampleDataset[1:500,],Labels<-Labels,"knn",n_genes_to_keep = 500)
tree<-Tree_Based_ML(ExampleDataset[1:500,],Labels=Labels,"C5.0",n_genes_to_keep = 500)

RankDataFrame=data.frame(matrix(nrow=500,ncol=4))
RankDataFrame$X1=scPNMF$Important_Features$InfoGene
RankDataFrame$X2=filtered_data$Important_Features
RankDataFrame$X3=Wrapper$Important_Features
RankDataFrame$X4=tree$Important_Features


match(unlist(RankDataFrame[,1]),rownames(ExampleDataset)
