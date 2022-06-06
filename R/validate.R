#' Run Silhouette
#'
#' This function calculates the average silhouette width of the clustering assignments made on the datasets provided
#'
#' @param dataList a data list of the data with clustering assignments
#' @param reduction_choosen The dimensionality reduction type of the data to find the calculated scores
#' @param method The clustering method used
#' @export
run_silhouette <- function(dataList, reduction_choosen = 'pca', method = 'kmeans')
{
  list_return = NULL
  for (i in (1:length(names(dataList))))
  {
    if(method == "kmeans"){
     c = paste("kmeans_cluster_", reduction_choosen, sep ="")
      x = dataList[[i]][[c]]
      y = Seurat::Embeddings(dataList[[i]], reduction = reduction_choosen)
      z <- cluster::silhouette(x[,1],dist(y, "euclidean"))
      x = summary(z)
      list_return = append(list_return, x$avg.width)
    }
    else if(method == "seurat")
    {
      x = dataList[[i]]@meta.data$seurat_clusters
      y = Seurat::Embeddings(dataList[[i]], reduction = reduction_choosen)
      z <- cluster::silhouette(x[,1],dist(y, "euclidean"))
      x = summary(z)
      list_return = append(list_return, x$avg.width)
    }
    else{
      stop("Invalid clustering method requested")
    }
  }
  return(list_return)
}

#' Run Dunn
#'
#' This function calculates the dunn index of the clustering assignments made on the datasets provided
#'
#' @param dataList a data list of the data with clustering assignments
#' @param reduction_choosen The dimensionality reduction type of the data to find the calculated scores
#' @param method The clustering method used
#' @export
run_dunn <- function(dataList, reduction_choosen,  method = 'kmeans'){

  list_return = NULL
  for (i in (1:length(names(dataList))))
  {
    if(method == 'kmeans'){
    c = paste("kmeans_cluster_", reduction_choosen, sep ="")
    y = dataList[[i]][[c]]
    x = Seurat::Embeddings(dataList[[i]], reduction = reduction_choosen)
    dunn = clValid::dunn(clusters = y[,1], Data = x)
    list_return = append(list_return, dunn)
    print(dunn)
    }
    else if(method == "seurat")
    {
      y = dataList[[i]]@meta.data$seurat_clusters
      x = Seurat::Embeddings(dataList[[i]], reduction = reduction_choosen)
      dunn = clValid::dunn(clusters = y[,1], Data = x)
      list_return = append(list_return, dunn)
      print(dunn)
    }
    else{
      stop("Invalid clustering method requested")
    }
    return(list_return)
  }
}


