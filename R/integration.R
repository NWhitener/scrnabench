#' Run CCA
#'
#' This function runs the data integration protocol detailed in the Seurat "Introduction to scRNA-seq
#' integration" found at https://satijalab.org/seurat/articles/integration_introduction.html. This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param dataList A data list of data sets to integrate using the cca protocol
#' @return data.combined A data list of the combined data from the cca protocol
#' @export
run_cca <- function(idx,seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(idx)
  dataList <- preprocess(dataList)
  #dataList <- run_gficf(dataList)
  dataList <- run_log(dataList)

  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- select_hvg(x)
  })

  features <- Seurat::SelectIntegrationFeatures(object.list = dataList)
  k.filter <- min(200, min(sapply(dataList, ncol)))
  data.anchors <- Seurat::FindIntegrationAnchors(object.list = dataList, anchor.features = features, k.filter=k.filter)
  data.combined <- Seurat::IntegrateData(anchorset = data.anchors)
  SeuratObject::DefaultAssay(data.combined) <- "integrated"

  data.combined <- scale_data(data.combined)
  data.combined <- run_pca(data.combined)
  data.combined <- run_umap(data.combnined)
  data.combined <- run_cluster(data.combined)

  return(data.combined)
}

#' Run Harmony
#'
#' This function runs the data integration protocol detailed in the Seurat ## FIND VING. AND LINK . This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param dataList A data list of data sets to integrate using the harmony protocol
#' @return data.combined A data list of the combined data from the harmony protocol
#' @export
run_harmony <- function(dataList, batch_name='ID', seed = 1)
{

  set.seed(seed)
  for (k in (1:length(names(dataList)))){
    print(k)
    dataList[[k]] <- harmony::RunHarmony(object = dataList[[k]], group.by.vars = batch_name)

  }
  return(dataList)
}

#' Run fastmnn
#'
#' This function runs the data integration protocol detailed in the Seurat ## FIND VING. AND LINK . This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param dataList A data list of data sets to integrate using the fastmnn protocol
#' @return data.combined A data list of the combined data from the fastmnn protocol
#' @export
run_fastmnn <- function(idx, batch_name = "ID", seed =1)
{

  ##REWRITE
  set.seed(seed)
  dataList <- extract_datasets(idx)
  dataList <- extract_common_genes(dataList)
  dataList <- merge_datasets(dataList)
  dataList <- preprocess(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_log(dataList) #LOG
  #dataList<- run_gficf(dataList)
  dataList <- select_hvg(dataList)
  dataList[[1]] <- SeuratWrapper::RunFastMNN(object.list = Seurat::SplitObject(dataList[[1]], split.by = batch_name))
  dataList <- run_umap(dataList)
  dataList <- run_cluster(dataList)
  SeuratObject::DefaultAssay(dataList) <- "mnn.reconstructed"
  return(dataList)

}


#' Run sctransform
#'
#' This function runs the data integration protocol detailed in the Seurat "Using sctransform in Seurtat"
#' found at https://satijalab.org/seurat/articles/sctransform_vignette.html . This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param dataList A data list of data sets to integrate using the sctransform protocol
#' @return data.combined A data list of the combined data from the sctransform protocol
#' @export
run_sctransform <- function(dataList, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(idx)
  dataList <- preprocess(dataList)
  dataList <- lapply(X = dataList, FUN = Seurat::SCTransform)
  features <- Seurat::SelectIntegrationFeatures(object.list = dataList, nfeatures=2000)
  dataList <- Seurat::PrepSCTIntegration(object.list = dataList, anchor.features = features)
  k.filter <- min(200, min(sapply(dataList, ncol)))
  data.anchors <- Seurat::FindIntegrationAnchors(object.list = dataList, normalization.method = "SCT", anchor.features = features, k.filter = k.filter)
  data.combined <- Seurat::IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
  SeuratObject::DefaultAssay(data.combined) <- "integrated"
  data.combined <- scale_data(data.combined)
  data.combined <- run_pca(data.combined)
  data.combined <- run_umap(data.combined)
  data.combined <- run_cluster(data.combined)
  return(data.combined)

}


#' Runs the Seurat Pipeline
#'
#' This function completes the Seurat Pipeline on a list of data. This process starts with finding the variable features (genes), then
#' scaling the data, running PCA, UMPA, and then find clusters.  This step can be done to test your benchamrk data sets
#' or for down stream analysis. For more information on the Seurat Pipeline
#' please see https://satijalab.org/seurat/index.html
#'
#' @param dataList A data list of Data sets in the Gene in Row and Cell in Columns format
#' @export
run_seurat <- function(dataList, seed =1)
{
  set.seed(seed)
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- Seurat::RunPCA(x, npcs = 30, verbose = FALSE)
    x <- Seurat::RunUMAP(x, reduction = "pca", dims = 1:30)
    x <- Seurat::FindNeighbors(x, reduction = "pca", dims = 1:30)
    x <- Seurat::FindClusters(x, resolution = 0.5)
  })
  return(dataList)
}
