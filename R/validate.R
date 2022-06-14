#' Run Silhouette
#'
#' This function calculates the average silhouette width of the clustering assignments made on the datasets provided
#'
#' @param dataList a data list of the data with clustering assignments
#' @param reductionChoosen The dimensionality reduction type of the data to find the calculated scores
#' @param method The clustering method used
#' @export
run_silhouette <- function(dataList, reductionChoosen = 'pca', method = 'kmeans')
{
  result = NULL
  for (i in (1:length(names(dataList))))
  {
    if(method == "kmeans"){
      reductionType <- paste("kmeans_cluster_", reductionChoosen, sep ="")
      clusters <- dataList[[i]][[reductionType]][,1]

    }
    else if(method == "seurat")
    {
      reductionType = paste('seurat_clusters')
      clusters <- as.numeric(dataList[[i]][[reductionType]]$seurat_clusters)
    }
    else{
      stop("Invalid clustering method requested")
    }
    if(length((unique(clusters))) < 2)
         {
          warningMessage = paste("Dataset ", datasets[i], " has only one cluster. Consider changing the clustering parameters", sep ='')
          warning(warningMessage)
          datasetResult = cbind(dataList[[i]]@meta.data$ID[1], NA)
          }
      else {
        cellEmbeddings <- Seurat::Embeddings(dataList[[i]], reduction = reductionChoosen)
        silhouetteScores <- summary(cluster::silhouette(clusters, dist(cellEmbeddings, "euclidean")))
        datasetResult = cbind(dataList[[i]]@meta.data$ID[1],  silhouetteScores$avg.width)
      }
    result = rbind(result, datasetResult)
    }
  return(result)
}

#' Run Dunn
#'
#' This function calculates the dunn index of the clustering assignments made on the datasets provided
#'
#' @param dataList a data list of the data with clustering assignments
#' @param reductionChoosen The dimensionality reduction type of the data to find the calculated scores
#' @param method The clustering method used
#' @export
run_dunn <- function(dataList, reductionChoosen,  method = 'kmeans'){
  result = NULL
  for (i in (1:length(names(dataList))))
  {
    if(method == "kmeans"){
      reductionType <- paste("kmeans_cluster_", reductionChoosen, sep ="")
      clusters <- dataList[[i]][[reductionType]][,1]

    }
    else if(method == "seurat")
    {
      reductionType = paste('seurat_clusters')
      clusters <- as.numeric(dataList[[i]][[reductionType]]$seurat_clusters)
    }
    else{
      stop("Invalid clustering method requested")
    }
    if(length((unique(clusters))) < 2)
    {
      warningMessage = paste("Dataset ", datasets[i], " has only one cluster. Consider changing the clustering parameters", sep ='')
      warning(warningMessage)
      datasetResult = cbind(dataList[[i]]@meta.data$ID[1], NA)
    }
    else {
      cellEmbeddings <- Seurat::Embeddings(dataList[[i]], reduction = reductionChoosen)
      dunn = clValid::dunn(clusters = clusters, Data = cellEmbeddings)
      datasetResult = cbind(dataList[[i]]@meta.data$ID[1],  dunn)
    }
    result = rbind(result, datasetResult)
  }
  return(result)
}


