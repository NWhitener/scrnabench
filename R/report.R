#' Create Internal Cluster Validation Report
#'
#' This function creates an internal validation report for clustering methods. This report includes the
#' Average Silhouette Width, as well as the Dunn Index as measures of Internal Validation. It also includes other useful information
#' such as the Number of Clusters and the Number of Singletons.
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
                          get_number_clusters(dataList,reductionType, method),
                          get_number_singletons(dataList,reductionType, method),
                          round(as.numeric(silhouetteScores[,2]), 2),
                          round(as.numeric(dunnScores[,2]), 2) )
    colNames <- c(colNames, paste('Number Clusters', toupper(reductionType), sep = ' '),paste("Number Singletons", toupper(reductionType), sep = ' '),
                  paste('Silhouette', toupper(reductionType), sep = ' '),
                  paste('Dunn', toupper(reductionType), sep = ' '))
  }

  resultsTable <- as.data.frame(resultsTable)
  colnames(resultsTable) = colNames
  return(resultsTable)
}

#' Create Metamorphic Testing Report
#'
#' This function creates a metamorphic testing report for clustering methods. This report includes which
#' reductions of the dataset passed the metamorphic test.
#'
#' @param reportList a list of internal clustering validation reports contain information from metamorphic testing
#' @return A report of the percent of reductions that passed the metamorphic test for each test and each dataset
#' @export
create_metamorphic_testing_report <- function(reportList)
{
  reportTable <- NULL
  if(is.list(dataList))
  {
    numberReductions <- sum(grepl('Number Clusters ', colnames(reportList[['original']])))
    if(numberReductions > 0)
    {
      reportList[['original']][is.na(reportList[['original']])] <- -10
      reportTable <- reportList[['original']][,1:2]
      numberItemsToCompare <- (ncol(reportList[['original']]) - 2 ) / numberReductions
      for (i in c(2:length(names(reportList))))
      {
        reportList[[i]][is.na(reportList[[i]])] <- -10
        matches <- (reportList[['original']] == reportList[[i]])
        testOutcomes <- rep(0, nrow(reportList[['original']]))
        for (j in seq(3, ncol(matches), by = numberItemsToCompare))
        {
          testOutcomes <- testOutcomes + ifelse(rowSums(matches[,j:(j + numberItemsToCompare - 1), drop =F]) == numberItemsToCompare, 1, 0)
        }
        reportTable[[names(reportList)[i]]] <- testOutcomes
      }
      reportTable[,3:ncol(reportTable)] = round(reportTable[,3:ncol(reportTable)] / numberReductions, 2)
    }
    else
    {
      stop("At least one reduction is required to create a report of metamorphic tests")
    }
  }
  else
  {
    stop("A list of datasets is required to create a report of metamorphic tests")
  }
  return(reportTable)
}


