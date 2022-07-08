#' Run CCA
#'
#' This function runs the data integration protocol detailed in the Seurat "Introduction to scRNA-seq
#' integration" found at https://satijalab.org/seurat/articles/integration_introduction.html. This function completes the cca process on
#' the datasets. Use run_cca_integration_workflow() for the entire workflow process
#'
#' @param dataList A list of data sets to integrate using the cca protocol
#' @return A data list of the combined data from the cca protocol
#' @export
run_cca <- function(dataList)
{
  if(is.list(dataList))
  {
    features <- Seurat::SelectIntegrationFeatures(object.list = dataList)
    kFilter <- min(200, min(sapply(dataList, ncol)))
    dataAnchors <- Seurat::FindIntegrationAnchors(object.list = dataList, anchor.features = features, k.filter=kFilter)
    print("Anchors Found")
    dataCombined <- list(Seurat::IntegrateData(anchorset = dataAnchors, k.weight= 50))
    names(dataCombined) <- c("Integrated")
  }
  else
  {
    stop("A data list of datasets is required to  run cca on the datasets")
  }
  return(dataCombined)
}

#' Run Harmony
#'
#' This function runs the data integration protocol found at https://github.com/satijalab/seurat-wrappers/blob/master/docs/harmony.md
#' This function completes the harmonization process on
#' the data sets. Use run_harmony_integration_workflow() for the entire workflow process
#'
#' @param dataList A list of data sets to integrate using the harmony protocol
#' @param batchName The annotation name to group data sets by, defaults to "ID"
#' @return  A list of the integrated data from the harmony protocol
#' @export
run_harmony <- function(dataList, batchName = 'ID')
{
  if(is.list(dataList))
  {
  for (i in (1:length(names(dataList)))){
    dataList[[i]] <- harmony::RunHarmony(object = dataList[[i]], group.by.vars = batchName)
  }
  }
  else
  {
    stop("A data list of datasets is required to run harmony on the datasets")
  }
  return(dataList)
}

#' Run fastmnn
#'
#' This function runs the data integration protocol detailed in the Seurat Wrappers protocol https://github.com/satijalab/seurat-wrappers/blob/master/docs/fast_mnn.md
#' This function completes the fastmnn process on
#' the datasets. Use run_fastmnn_integration_workflow() for the entire workflow process
#'
#' @param dataList A  list of data sets to integrate using the fastmnn protocol
#' @param batchName The annotation name to group data sets by, defaults to "ID"
#' @return A data list of the combined data from the fastmnn protocol
#' @export
run_fastmnn <- function(dataList, batchName = "ID")
{
  if(is.list(dataList))
  {
  for (i in 1:length(names(dataList))){
    dataList[[i]] <- SeuratWrappers::RunFastMNN(object.list = Seurat::SplitObject(dataList[[i]], split.by = batchName))
  }
  }
  else
  {
    stop("A data list of datasets is required to run fastmnn on the datasets")
  }
  return(dataList)

}


#' Run sctransform
#'
#' This function runs the data integration protocol detailed in the Seurat "Using sctransform in Seurtat"
#' found at https://satijalab.org/seurat/articles/sctransform_vignette.html . This function completes the sctransform process on
#' the datasets. Use run_sctransform_workflow() for the entire workflow process
#'
#' @param dataList A list of data sets to integrate using the sctransform protocol
#' @param numFeatures The number of Features to use, defaults to 2000
#' @return  A data list of the combined data from the sctransform protocol
#' @export
run_sctransform <- function(dataList, numFeatures = 2000)
{

  if(is.list(dataList))
  {
  dataList <- lapply(X = dataList, FUN = Seurat::SCTransform)
  features <- Seurat::SelectIntegrationFeatures(object.list = dataList, nfeatures=numFeatures)
  dataList <- Seurat::PrepSCTIntegration(object.list = dataList, anchor.features = features)
  kFilter <- min(200, min(sapply(dataList, ncol)))
  dataAnchors <- Seurat::FindIntegrationAnchors(object.list = dataList, normalization.method = "SCT", anchor.features = features, k.filter = kFilter)
  dataCombined <- list(Seurat::IntegrateData(anchorset = dataAnchors, normalization.method = "SCT"))
  names(dataCombined)<-c("Integrated")
  }
  else
  {
    stop("A data list of datasets is required to run sctransform on the datasets")
  }
  return(dataCombined)
}



