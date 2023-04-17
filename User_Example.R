library(SeuratData)
library(scrnabench)
library(Seurat)

filter_data_real <- function(dataList)
{
  nCount_RNA = NULL
  nFeature_RNA = NULL
  percent.mt = NULL
  if(is.list(dataList))
  {
    for (i in (1:length(dataList)))
    {
      dataList[[i]][["percent.mt"]] <- Seurat::PercentageFeatureSet(dataList[[i]], pattern = "^MT-")
      Total_mRNAs <- dataList[[i]][["nCount_RNA"]]$nCount_RNA
      mupper_bound <- 10^(mean(log10(Total_mRNAs)) + 2*stats::sd(log10(Total_mRNAs)))
      mlower_bound <- 10^(mean(log10(Total_mRNAs)) - 2*stats::sd(log10(Total_mRNAs)))
      Total_Genes <- dataList[[i]][["nFeature_RNA"]]$nFeature_RNA
      gupper_bound <- 10^(mean(log10(Total_Genes)) + 2*stats::sd(log10(Total_Genes)))
      glower_bound <- 10^(mean(log10(Total_Genes)) - 2*stats::sd(log10(Total_Genes)))
      dataList[[i]] <- subset(x = dataList[[i]], subset = nFeature_RNA > glower_bound & nFeature_RNA < gupper_bound &
                                nCount_RNA > mlower_bound & nCount_RNA < mupper_bound & percent.mt < 10)
    }
  }
  else{
    stop("A data list of datasets is required to preprocess datasets")
  }
  
  return(dataList)
}

run_harmony_integration_workflow_real <- function(dataList, batchName)
{
  if(is.list(dataList))
  {
    dataList <- filter_data_real(dataList)
    dataList <- run_log(dataList)
    dataList <- select_hvg(dataList)
    dataList <- scale_data(dataList)
    dataList <- run_pca(dataList, numComponents = 10)
    dataList <- run_harmony(dataList, batchName = batchName)
    dataList <- run_umap(dataList, reductionType = 'harmony', numDimensions = 10)
  }
  else
  {
    stop("A data list of datasets is required to apply the harmony workflow to the datasets")
  }
  return(dataList)
}

run_fastmnn_integration_workflow_real <- function(dataList, batchName)
{
  if(is.list(dataList))
  {
    dataList <- filter_data_real(dataList)
    dataList <- run_log(dataList)
    dataList <- select_hvg(dataList)
    dataList <- run_fastmnn(dataList, batchName = batchName)
    dataList <- run_umap(dataList, reductionType = 'mnn', numDimensions = 10)
  }
  return(dataList)
}

run_cca_integration_workflow_real <- function(dataList, batchName)
{
  if(is.list(dataList))
  {
    dataList <- SplitObject(dataList$Integrated, split.by = batchName)
    dataList <- filter_data_real(dataList)
    dataList <- run_log(dataList)
    dataList <- select_hvg(dataList)
    dataList <- run_cca(dataList)
    dataList <- scale_data(dataList)
    dataList <- run_pca(dataList, numComponents = 10)
    dataList <- run_umap(dataList, reductionType = 'pca', numDimensions = 10)
  }
  else
  {
    stop("A data list of datasets is required to apply the cca workflow to the datasets")
  }
  return(dataList)
}

run_sctransform_integration_workflow_real <- function(dataList, batchName)
{
  if(is.list(dataList))
  {
    dataList <- SplitObject(dataList$Integrated, split.by = batchName)
    dataList <- filter_data_real(dataList)
    dataList <- run_sctransform(dataList)
    dataList <- scale_data(dataList)
    
    dataList <- run_pca(dataList, numComponents = 10)
    dataList <- run_umap(dataList, reductionType = 'pca', numDimensions = 10)
  }
  else
  {
    stop("A data list of datasets is required to apply the cca workflow to the datasets")
  }
  return(dataList)
}

directory = "../data/"
dataList <- c("panc8", "pbmcsca", "ifnb")
for (i in dataList)
{
  data <- LoadData(i)
  data <- UpdateSeuratObject(object = data)
  data <- list(data)
  names(data) <- c('Integrated')

  subdirectory = paste(directory, 'real_data/', i, sep="")
  if(!file.exists(subdirectory))
  {
    dir.create(subdirectory)
  }
  
  batchId <- NULL  
  if(i == "panc8")
  {
    data$Integrated <- subset(x = data$Integrated, subset = tech == "fluidigmc1", invert = TRUE)
    batch <- "tech"
  }
  else if(i == "pbmcsca")
  {
    data$Integrated <- subset(x = data$Integrated, subset = Method == "CEL-Seq2", invert = TRUE)
    data$Integrated <- subset(x = data$Integrated, subset = Method == "Smart-seq2", invert = TRUE)
    batch <- "Method"
  }
  else if(i == "ifnb")
  {
    batch <- "stim"
  }
  
  integratedData <- run_harmony_integration_workflow_real(data, batchName = batch)
  umapEmbeddings <- Seurat::Embeddings(integratedData$Integrated, reduction = "umap")
  write.csv(umapEmbeddings, file=paste(subdirectory, "/harmony_umap.csv", sep=""))
  write.csv(integratedData$Integrated@meta.data, file=paste(subdirectory, "/harmony_labels.csv", sep=""))

  integratedData <- run_fastmnn_integration_workflow_real(data, batchName = batch)
  umapEmbeddings <- Seurat::Embeddings(integratedData$Integrated, reduction = "umap")
  write.csv(umapEmbeddings, file=paste(subdirectory, "/fastmnn_umap.csv", sep=""))
  write.csv(integratedData$Integrated@meta.data, file=paste(subdirectory, "/fastmnn_labels.csv", sep=""))
  
  integratedData <- run_cca_integration_workflow_real(data, batchName = batch)
  umapEmbeddings <- Seurat::Embeddings(integratedData$Integrated, reduction = "umap")
  write.csv(umapEmbeddings, file=paste(subdirectory, "/cca_umap.csv", sep=""))
  write.csv(integratedData$Integrated@meta.data, file=paste(subdirectory, "/cca_labels.csv", sep=""))
  
  integratedData <- run_sctransform_integration_workflow_real(data, batchName = batch)
  umapEmbeddings <- Seurat::Embeddings(integratedData$Integrated, reduction = "umap")
  write.csv(umapEmbeddings, file=paste(subdirectory, "/sctransform_umap.csv", sep=""))
  write.csv(integratedData$Integrated@meta.data, file=paste(subdirectory, "/sctransform_labels.csv", sep=""))
  
}




