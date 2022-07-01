#' Create Internal Cluster Validation Report
#'
#' This function creates an internal validation report for clustering methods. This report includes the
#' Average Silhouette Width, as well as the Dunn Index as measures of Internal Validation.
#'
#' @param dataList A list of datasets with clustering completed
#' @param method The type of clustering method to validate, defaults to kmeans
#' @return A Internal Cluster Validation Report table with the Number of clusters, Average Silhouette Width,
#' and Dunn Index as measures of validation
#' @export
create_internal_cluster_validation_report <- function(dataList, method = 'kmeans')
{
  resultsTable <- NULL
  colNames <- c('Dataset', 'Number of Cells')

  for(i in 1:length(names(dataList)))
  {
    resultsTable <- rbind(resultsTable, cbind(names(dataList[i]), ncol(dataList[[i]])))
  }

  reductions <- names(dataList[[1]]@reductions)

  for (reductionType in reductions)
  {
    annotationField <- toupper(paste(method, '_cluster_', reductionType, sep = ''))
    silhouetteScores <- run_silhouette(dataList, reductionType = reductionType, method = method)
    dunnScores <- run_dunn(dataList, reductionType = reductionType, method = method)
    resultsTable <- cbind(resultsTable,
                          get_num_clusters(dataList, reductionType, method),
                          round(as.numeric(silhouetteScores[,2]), 2),
                          round(as.numeric(dunnScores[,2]), 2))
    colNames <- c(colNames, paste('Number Clusters', toupper(reductionType), sep = ' '),
                  paste('Silhouette', toupper(reductionType), sep = ' '),
                  paste('Dunn', toupper(reductionType), sep = ' '))
    }

  resultsTable <- as.data.frame(resultsTable)
  colnames(resultsTable) = colNames
  return(resultsTable)
}



get_num_clusters <- function(dataList, reductionType= 'pca', method = 'kmeans')
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
