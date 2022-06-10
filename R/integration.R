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
    features <- SelectIntegrationFeatures(object.list = dataList)
    k.filter <- min(200, min(sapply(dataList, ncol)))
    data.anchors <- FindIntegrationAnchors(object.list = dataList, anchor.features = features, k.filter=k.filter)
    data.combined <- IntegrateData(anchorset = data.anchors)
    DefaultAssay(data.combined) <- "integrated"

    return(list(data.combined))
}

#' Run Harmony
#'
#' This function runs the data integration protocol detailed in the Seurat ## FIND VING. AND LINK . This function completes the harmonization process on
#' the datasets. Use run_harmony_workflow() for the entire workflow process
#'
#' @param dataList A data list of data sets to integrate using the harmony protocol
#' @return data.combined A data list of the combined data from the harmony protocol
#' @export
run_harmony <- function(dataList, batch_name='ID')
{


  for (k in (1:length(names(dataList)))){
    print(k)
    dataList[[k]] <- harmony::RunHarmony(object = dataList[[k]], group.by.vars = batch_name)

  }
  return(dataList)
}

#' Run fastmnn
#'
#' This function runs the data integration protocol detailed in the Seurat ## FIND VING. AND LINK . This function completes the fastmnn process on
#' the datasets. Use run_fastmnn_workflow() for the entire workflow process
#'
#' @param dataList A data list of data sets to integrate using the fastmnn protocol
#' @return data.combined A data list of the combined data from the fastmnn protocol
#' @export
run_fastmnn <- function(idx, batch_name = "ID", seed =1)
{
  for (k in (1:length(names(dataList)))){
    print(k)
    dataList[[k]] <- SeuratWrapper::RunFastMNN(object.list = Seurat::SplitObject(dataList[[k]], split.by = batch_name))

  }
  return(dataList)

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
  dataList <- lapply(X = dataList, FUN = SCTransform)
  features <- SelectIntegrationFeatures(object.list = dataList, nfeatures=2000)
  dataList <- PrepSCTIntegration(object.list = dataList, anchor.features = features)
  k.filter <- min(200, min(sapply(dataList, ncol)))
  data.anchors <- FindIntegrationAnchors(object.list = dataList, normalization.method = "SCT", anchor.features = features, k.filter = k.filter)
  data.combined <- IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
  DefaultAssay(data.combined) <- "integrated"

  return(list(data.combined))
}



