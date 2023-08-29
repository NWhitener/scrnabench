#' Run Silhouette
#'
#' This function calculates the average silhouette width of the clustering assignments made on the data sets provided.
#' Uses the clusterCrit package implementation.
#'
#'
#' @param dataList a data list of the data with clustering assignments
#' @param reductionType The dimensionality reduction type of the data to find the calculated scores,
#' defaults to PCA
#' @param method The clustering method used, defaults to Kmeans
#' @param sampling The size of the data sample that should be used
#' @param t The iteration of resampling the dataset and re-calculate silhouette scores
#' @return A table of the Silhouette Scores, with the data set name and the score
#' @export
run_silhouette <- function(dataList, method = 'kmeans', reductionType = 'pca', sampling = 8000, t = 2)
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
      if(length(cellEmbeddings[,1]) > sampling)
      {
        silhouetteScore = 0
        for(i in (1:t))
        {
          sampleIndices <- sample(nrow(cellEmbeddings), sampling)
          samplesData <- cellEmbeddings[sampleIndices, ]
          samplesClusters <- clusters[sampleIndices]
          sampleSil <- mean(cluster::silhouette(as.integer(samplesClusters), dist(samplesData))[,3])
          silhouetteScore <-silhouetteScore +  sampleSil
        }
        silhouetteScore <- silhouetteScore / t
      }
      else
      {
        silhouetteScore <- mean(cluster::silhouette(as.integer(clusters), dist(cellEmbeddings))[,3])
      }
      resultsTable = rbind(resultsTable, c(names(dataList[i]),  silhouetteScore))
    }
  }
  return(resultsTable)
}

#' Run Dunn Index
#'
#' This function calculates the Dunn index of the clustering assignments made on the datasets provided. Uses the
#' clusterCrit implementation of the Dunn Index
#'
#' @param dataList a data list of the data with clustering assignments
#' @param reductionType The dimensionality reduction type of the data to find the calculated scores,
#' defaults to PCA
#' @param method The clustering method used, defaults to Kmeans
#' @param sampling The size of the data sample that should be used
#' @param t The iteration of resampling the dataset and re-calculate silhouette scores
#' @return A table of the Dunn Index Scores, with the data set name and the score
#' @export
run_dunn <- function(dataList, reductionType = 'pca',  method = 'kmeans', sampling = 8000, t = 2)
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
      if(length(cellEmbeddings[,1]) > sampling)
      {
        dunnScore = 0
        for(i in (1:t))
        {
          sampleIndices <- sample(nrow(cellEmbeddings), sampling)
          samplesData <- cellEmbeddings[sampleIndices, ]
          samplesClusters <- clusters[sampleIndices]
          sampleDunn <- clValid::dunn(Data = samplesData, clusters = as.integer(samplesClusters))
          dunnScore <- dunnScore +  sampleDunn
        }
        dunnScore <- dunnScore / t
      }
      else
      {
        dunnScore =  clValid::dunn(Data = cellEmbeddings, clusters = as.integer(clusters))
      }
      
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



