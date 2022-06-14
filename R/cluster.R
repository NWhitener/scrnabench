#' Run Kmeans
#'
#' This function runs Kmeans Clustering. It uses the simple Kmeans algorithm with
#' a max number of iterations set to be 100.
#'
#' @param dataList A list of datasets to apply the Kmeans algorithm too
#' @param k The number of cluster centers to use, defaults to 10
#' @return The dataList with kmeans clusters stored in the meta data
#' @export
run_kmeans <- function(dataList, k=10, reductionChoosen = 'pca')
{
  for (i in (1:length(names(dataList))))
  {
    clustData = Seurat::Embeddings(dataList[[i]], reduction = reductionChoosen)
    meta = dataList[[i]]@meta.data
    meta[paste("kmeans_cluster_", reductionChoosen, sep = "")]<- stats::kmeans(clustData,  k, iter.max = 100)$cluster
    print(i)
    metaFix <- subset(meta, select = c(paste("kmeans_cluster_", reductionChoosen, sep = "")))
    dataList[[i]] <- Seurat::AddMetaData(dataList[[i]], metaFix)

  }

  return(dataList)
}

#' Complete Clustering Steps
#'
#' This functions runs the Seurat FindNeighbors and FindClusters function on a list of data sets.
#' This functions assumes that genes are in rows and
#'  cells are in columns. The FindNeighbors reduction is set to "pca" by default and uses the first 30 dimensions
#'  The FindClusters resolution is set to 0.5
#'
#' @param dataList A data list of data sets
#' @return A data list with clustering completed completed on the features
#' @export
run_seurat_cluster <- function(dataList, reductionChoosen = "pca", resolutionGiven = 0.5)
{
  for (i in (1:length(names(dataList)))){
    dataList[[i]] <- Seurat::FindNeighbors(dataList[[i]], reduction = reductionChoosen, dims = 1:30)
    dataList[[i]]<- Seurat::FindClusters(dataList[[i]], resolution = resolutionGiven)
  }
  return(dataList)
}
