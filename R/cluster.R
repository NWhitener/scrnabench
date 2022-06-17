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

  ##REFACTOR TO THE ADD ANNOTATION FUNCTION

  if(is.list(dataList))
  {
  for (i in (1:length(names(dataList))))
  {
    clustData = Seurat::Embeddings(dataList[[i]], reduction = reductionChoosen)
    meta = dataList[[i]]@meta.data
    meta[paste("kmeans_cluster_", reductionChoosen, sep = "")] <- stats::kmeans(clustData,  k, iter.max = 100)$cluster
    print(i)
    metaFix <- subset(meta, select = c(paste("kmeans_cluster_", reductionChoosen, sep = "")))
    dataList[[i]] <- Seurat::AddMetaData(dataList[[i]], metaFix)

  }
    return(dataList)
  }
  else
    {
      stop("A data list of datasets is required to kmeans cluster datasets")
    }

}

#' Complete Clustering Steps
#'
#' This functions runs the Seurat FindNeighbors and FindClusters function on a list of data sets.
#' This functions assumes that genes are in rows and
#'  cells are in columns. The FindNeighbors reduction is set to "pca" by default and uses the first 30 dimensions
#'  The FindClusters resolution is set to 0.5
#'
#' @param dataList A data list of data sets
#' @param reductionChoosen The dimensionality reduction type that is wanted
#' @param resolutionGiven the resolution of the seurat clustering method
#' @param numComponents the number of components to use
#' @return A data list with clustering completed completed on the features
#' @export
run_seurat_cluster <- function(dataList, reductionChoosen = "pca", resolutionGiven = 0.5, numComponents = 10)
{

## rewrite, go through objects find clusters
  if(is.list(dataList))
  {
  for (i in (1:length(names(dataList)))){
    dataList[[i]] <- Seurat::FindNeighbors(dataList[[i]], reduction = reductionChoosen, dims=1:numComponents)
    dataList[[i]]<- Seurat::FindClusters(dataList[[i]], resolution = resolutionGiven)
    ##Find clusters
    ##store them into list of lists

  }

## call add annotation
  return(dataList)


    }
  else
  {
    stop("A data list of datasets is required to use seurat clustering on the datasets")
  }
}
#' Find the number of clusters
#'
#' This functions finds the number of clusters in the Seurat Clustering method, or through the kmeans clustering
#' method.
#'
#' @param dataList A data list of data sets
#' @param reductionChoosen The dimensionality reduction type that is wanted
#' @param clusteringType the type of clustering you would like to find the number of
#' @return A results list with the number of clusters
#' @export
find_num_cluster <- function(dataList, clusteringType = 'kmeans', reductionType= 'pca')
{
  if(is.list(dataList))
  {
    result = NULL
    for (i in (1:length(names(dataList)))){
    if(clusteringType == 'kmeans')
    {
      reductionType <- paste("kmeans_cluster_", reductionChoosen, sep ="")
      clusters <- dataList[[i]][[reductionType]][,1]
    }
    else if(clusteringType == 'seurat')
    {
      reductionType = paste('seurat_clusters')
      clusters <- as.numeric(dataList[[i]][[reductionType]]$seurat_clusters)
    }
    numClust = length(unique(clusters))
    result = append(result, numClust)
    }
    return(result)


  }
}
