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
    temp_list= NULL
    temp_list = cbind(temp_list, dataList[[i]]@meta.data$ID[1])
    if(method == "kmeans"){
     c = paste("kmeans_cluster_", reduction_choosen, sep ="")
      x = dataList[[i]][[c]]
      y = Seurat::Embeddings(dataList[[i]], reduction = reduction_choosen)
      z <- cluster::silhouette(x[,1],dist(y, "euclidean"))
      x = summary(z)
      temp_list = cbind(temp_list,  x$avg.width)
    }
    else if(method == "seurat")
    {
      c = paste('seurat_clusters')
      x = dataList[[i]][[c]]
      x = as.numeric(x$seurat_clusters)
      y = Seurat::Embeddings(dataList[[i]], reduction = reduction_choosen)
      z <- cluster::silhouette(x, dist(y, "euclidean"))
      q = summary(z)
      print(i)
      print(class(q))
      temp_list = cbind(temp_list, q$avg.width)
    }
    else{
      stop("Invalid clustering method requested")
    }
    list_return = rbind(list_return, temp_list)
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
    temp_list= NULL
    temp_list = cbind(temp_list, dataList[[i]]@meta.data$ID[1])
    if(method == 'kmeans'){
    c = paste("kmeans_cluster_", reduction_choosen, sep ="")
    y = dataList[[i]][[c]]
    x = Seurat::Embeddings(dataList[[i]], reduction = reduction_choosen)
    dunn = clValid::dunn(clusters = y[,1], Data = x)
    temp_list = cbind(temp_list, dunn)

    }
    else if(method == "seurat")
    {
      c = paste('seurat_clusters')
      y = dataList[[i]][[c]]
      y <- as.numeric(y$seurat_clusters)
      x = Seurat::Embeddings(dataList[[i]], reduction = reduction_choosen)
      dunn = clValid::dunn(clusters = y, Data = x)
      temp_list = cbind(temp_list, dunn)
    }
    else{
      stop("Invalid clustering method requested")
    }
    list_return = rbind(list_return, temp_list)
  }
  return(list_return)
}


