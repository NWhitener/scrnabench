#' Run Kmeans
#'
#' This function runs Kmeans Clustering. It uses the simple Kmeans algorithm with
#' a max number of iterations set to be 100 by default. The number of clusters is set by the users and defaults
#' to 10 clusters. This is the R package stats implementation of Kmeans
#'
#' @param dataList A list of data sets to be clustered by Kmeans
#' @param reductionType The type of dimensionly reduced data that should be used for clustering
#' defaults to PCA
#' @param numberClusters The number of cluster centers to use, defaults to 10
#' @param iterationsMax The max number of iterations that should be used in the kmeans algorithm
#' defaults to 100
#' @return The dataList with kmeans clusters stored in the meta data
#' @export
run_kmeans_clustering <- function(dataList, reductionType = 'pca', numberClusters = 10, iterationsMax = 100)
{
  if(is.list(dataList))
  {
    for (i in (1:length(names(dataList))))
    {
      cellEmbeddings = Seurat::Embeddings(dataList[[i]], reduction = reductionType)
      annotationField <- toupper(paste('kmeans_cluster_', reductionType, sep=''))
      dataList[[i]][[annotationField]] <- stats::kmeans(cellEmbeddings, numberClusters, iter.max = iterationsMax)$cluster
    }
  }
  else
  {
    stop("A data list of datasets is required to run k-means clustering.")
  }

  return(dataList)
}

#' Seurat Clustering
#'
#' This functions runs the Seurat FindNeighbors and FindClusters function on a list of data sets.
#' This functions assumes that genes are in rows and
#' cells are in columns. The FindNeighbors reduction is set to "pca" by default and uses the first 10 dimensions
#' The FindClusters resolution is set to 0.5 by default, but is controlable by the user
#'
#' @param dataList A list of data sets to be clustered
#' @param reductionType The type of dimensionly reduced data that should be used for clustering, defaults to PCA
#' @param resolutionValue The resolution for the Seurat clustering method, defaults to 0.5
#' @param numberComponents the number of components to use, defaults to 10
#' @return A data list with Seurat clustering completed
#' @export
run_seurat_clustering <- function(dataList, reductionType = 'pca', resolutionValue = 0.5, numberComponents = 10)
{

  if(is.list(dataList))
  {
    for (i in (1:length(names(dataList))))
    {
      dataList[[i]] <- Seurat::FindNeighbors(dataList[[i]], reduction = reductionType, dims = 1:numberComponents)
      dataList[[i]]<- Seurat::FindClusters(dataList[[i]], resolution = resolutionValue)
      annotationField <- toupper(paste('seurat_cluster_', reductionType, sep=''))
      dataList[[i]][[annotationField]] <- as.numeric(dataList[[i]][['seurat_clusters']][,1])
      dataList[[i]] <- Seurat::DietSeurat(dataList[[i]], data = T,
                                          counts = F, scale.data = F, features = Seurat::VariableFeatures(object = dataList[[i]]),
                                          dimreducs = names(dataList[[i]]@reductions), graphs = NULL)
    }
  }
  else
  {
    stop("A list of datasets is required to run seurat clustering.")
  }

  return(dataList)
}



#' Get the number of clusters
#'
#' This functions finds the number of clusters in the Seurat Clustering method, or through the kmeans clustering
#' method.
#'
#' @param dataList A data list of data sets
#' @param method The type of clustering, defaults to kmeans
#' @param reductionType The type of dimensionly reduced data that should be used for clustering,
#' defaults to PCA
#' @return A list with the number of clusters
#' @export
get_number_clusters <- function(dataList, reductionType= 'pca', method = 'kmeans')
{
  if(is.list(dataList))
  {
    numClusters = NULL
    for(i in 1:length(names(dataList)))
    {
      annotationField <- toupper(paste(method, '_cluster_', reductionType, sep = ''))
      numClusters = append(numClusters, length(unique(dataList[[i]][[annotationField]][,1])))
    }
    return(numClusters)
  }
  else
  {
    stop("A data list of datasets is required to get the number of clusters in the datasets")
  }

}


get_number_singletons <- function(dataList, reductionType = 'pca', method = 'kmeans')
{
  if(is.list(dataList))
  {
    numSingletons = NULL
    for(i in 1:length(names(dataList)))
    {
      annotationField <- toupper(paste(method, '_cluster_', reductionType, sep = ''))

      clusters <- dataList[[i]][[annotationField]][,1]
      singletonClusters <- which(data.frame(table(clusters))$Freq == 1)
      numSingletons = append(numSingletons, length(singletonClusters))
    }
    print(numSingletons)
    return(numSingletons)
  }
  else
  {
    stop("A data list of datasets is required to get the number of clusters in the datasets")
  }

}
