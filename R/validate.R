#' Run Silhouette
#'
#' This function calculates the average silhouette width of the clustering assignments made on the data sets provided
#' Uses the cluster package implementation.
#'
#'
#' @param dataList a data list of the data with clustering assignments
#' @param reductionType The dimensionality reduction type of the data to find the calculated scores,
#' defaults to PCA
#' @param method The clustering method used, defaults to Kmeans
#' @return A table of the Silhouette Scores, with the data set name and the score
#' @export
run_silhouette <- function(dataList, method = 'kmeans', reductionType = 'pca')
{
  resultsTable = NULL
  annotationField <- toupper(paste(method, '_cluster_', reductionType, sep=''))
  for (i in (1:length(names(dataList))))
  {
    clusters <- dataList[[i]][[annotationField]][,1]
    singletonClusters <- which(data.frame(table(clusters))$Freq == 1)
    if(length(unique(clusters)) < 2 | length(singletonClusters) > 0)
    {
      warningMessage = paste("Dataset ", names(dataList[i]), " has only one cluster or contains singletons. Consider changing the clustering parameters", sep ='')
      warning(warningMessage)
      resultsTable = rbind(resultsTable, c(names(dataList[i]), NA))
    }
    else
    {
      cellEmbeddings <- Seurat::Embeddings(dataList[[i]], reduction = reductionType)
      silhouetteScore <- clusterCrit::intCriteria(cellEmbeddings, as.integer(clusters), c("Silhouette"))
      resultsTable = rbind(resultsTable, c(names(dataList[i]),  silhouetteScore))
    }
  }
  return(resultsTable)
}

#' Run Dunn Index
#'
#' This function calculates the Dunn index of the clustering assignments made on the datasets provided. Uses the
#' clVaild implementation of the Dunn Index
#'
#' @param dataList a data list of the data with clustering assignments
#' @param reductionType The dimensionality reduction type of the data to find the calculated scores,
#' defaults to PCA
#' @param method The clustering method used, defaults to Kmeans
#' #' @return A table of the Dunn Index Scores, with the data set name and the score
#' @export
run_dunn <- function(dataList, reductionType = 'pca',  method = 'kmeans')
{
  resultsTable = NULL
  annotationField <- toupper(paste(method, '_cluster_', reductionType, sep=''))

  for (i in (1:length(names(dataList))))
  {
    clusters <- dataList[[i]][[annotationField]][,1]
    if(length(unique(clusters)) < 2)
    {
      warningMessage = paste("Dataset ", names(dataList[i]), " has only one cluster. Consider changing the clustering parameters", sep ='')
      warning(warningMessage)
      resultsTable = rbind(resultsTable, c(names(dataList[i]), NA))
    }
    else
    {
      cellEmbeddings <- Seurat::Embeddings(dataList[[i]], reduction = reductionType)
      dunnScore =  clusterCrit::intCriteria(cellEmbeddings, as.integer(clusters), c("Dunn"))
      resultsTable = rbind(resultsTable, c(names(dataList[i]),  dunnScore))
    }
  }
  return(resultsTable)
}

#' Run ARI
#'
#' Computes the Adjusted Rand Index of the cluster memberships. Uses the aricode implementation of ARI
#'
#' @param dataList A list of data sets with cluster membership information
#' @param groundTruths The true cell cluster memberships
#' @param reductionType The dimensionality reduction type of the data to find the calculated scores,
#' defaults to PCA
#' @param method The clustering method used, defaults to Kmeans
#' @return A table of the ARI Scores, with the data set name and the score
#' @export
run_ari <- function(dataList, groundTruths, reductionType = 'pca',  method = 'kmeans')
{
  resultsTable = NULL
  annotationField <- toupper(paste(method, '_cluster_', reductionType, sep=''))
  for (i in (1:length(names(dataList))))
  {
    clusters <- dataList[[i]][[annotationField]][,1]
    ariScore <- aricode::ARI(clusters, groundTruths[[i]])
    resultsTable = rbind(resultsTable, c(names(dataList[i]), ariScore))
  }
  return(resultsTable)
}






