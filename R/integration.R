#' Run CCA
#'
#' This function runs the data integration protocol detailed in the Seurat "Introduction to scRNA-seq
#' integration" found at https://satijalab.org/seurat/articles/integration_introduction.html. This function completes the cca process on
#' the datasets. Use run_cca_workflow() for the entire workflow process
#'
#' @param dataList A data list of data sets to integrate using the cca protocol
#' @return data.combined A data list of the combined data from the cca protocol
#' @export
run_cca <- function(dataList)
  {

  if(is.list(dataList))
  {
    features <- Seurat::SelectIntegrationFeatures(object.list = dataList)
    k.filter <- min(200, min(sapply(dataList, ncol)))
    data.anchors <- Seurat::FindIntegrationAnchors(object.list = dataList, anchor.features = features, k.filter=k.filter)
    data.combined <- Seurat::IntegrateData(anchorset = data.anchors)
    Seurat::DefaultAssay(data.combined) <- "integrated"
    merged<-list(data.combined)
    names(merged)<-c("Integrated")
    return(merged)
  }
  else
  {
    stop("A data list of datasets is required to  run cca on the datasets")
  }
  }

#' Run Harmony
#'
#' This function runs the data integration protocol detailed in the Seurat ## FIND VING. AND LINK . This function completes the harmonization process on
#' the datasets. Use run_harmony_workflow() for the entire workflow process
#'
#' @param dataList A data list of data sets to integrate using the harmony protocol
#' @return data.combined A data list of the combined data from the harmony protocol
#' @export
run_harmony <- function(dataList, batchName='ID')
{

  if(is.list(dataList))
  {
  for (k in (1:length(names(dataList)))){
    print(k)
    dataList[[k]] <- harmony::RunHarmony(object = dataList[[k]], group.by.vars = batchName)

  }
  return(dataList)
  }
  else
  {
    stop("A data list of datasets is required to run harmony on the datasets")
  }
}

#' Run fastmnn
#'
#' This function runs the data integration protocol detailed in the Seurat ## FIND VING. AND LINK . This function completes the fastmnn process on
#' the datasets. Use run_fastmnn_workflow() for the entire workflow process
#'
#' @param dataList A data list of data sets to integrate using the fastmnn protocol
#' @return data.combined A data list of the combined data from the fastmnn protocol
#' @export
run_fastmnn <- function(dataList, batchName = "ID", seed =1)
{

  if(is.list(dataList))
  {
  for (k in 1:length(names(dataList))){
    print(k)
    print(class(dataList[[k]]))
    dataList[[k]] <- SeuratWrappers::RunFastMNN(object.list = Seurat::SplitObject(dataList[[k]], split.by = batchName))

  }
  return(dataList)}
  else
  {
    stop("A data list of datasets is required to run fastmnn on the datasets")
  }


}


#' Run sctransform
#'
#' This function runs the data integration protocol detailed in the Seurat "Using sctransform in Seurtat"
#' found at https://satijalab.org/seurat/articles/sctransform_vignette.html . This function completes the sctransform process on
#' the datasets. Use run_sctransform_workflow() for the entire workflow process
#'
#' @param dataList A data list of data sets to integrate using the sctransform protocol
#' @return data.combined A data list of the combined data from the sctransform protocol
#' @export
run_sctransform <- function(dataList)
{

  if(is.list(dataList))
  {
  dataList <- lapply(X = dataList, FUN = Seurat::SCTransform)
  features <- Seurat::SelectIntegrationFeatures(object.list = dataList, nfeatures=2000)
  dataList <- Seurat::PrepSCTIntegration(object.list = dataList, anchor.features = features)
  k.filter <- min(200, min(sapply(dataList, ncol)))
  data.anchors <- Seurat::FindIntegrationAnchors(object.list = dataList, normalization.method = "SCT", anchor.features = features, k.filter = k.filter)
  data.combined <- Seurat::IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
  Seurat::DefaultAssay(data.combined) <- "integrated"
  merged<-list(data.combined)
  names(merged)<-c("Integrated")
  return(merged)
  }
  else
  {
    stop("A data list of datasets is required to run sctransform on the datasets")
  }


}



