#' Preturb Data sets
#'
#' This function applies a specific perturbation to the data-list based on the perturbation type. For a description of the
#' perturbations please see the individual functions
#'
#' @param dataList A list of data set that should be perturbed
#' @param perturbationType The number of the perturbation that is too be applied
#' @return a data list of the perturbed datasets
#' @export
perturb_datasets <- function(dataList, perturbationType = 1)
{
  dataList <- switch(
    perturbationType,
    '1'= permute_columns(dataList),
    '2'= modify_gene_counts(dataList),
    '3'= add_duplicate_cells(dataList),
    '4'= permute_rows(dataList),
    '5'= add_zero_variance_gene_counts(dataList),
    '6'= flip_gene_counts(dataList)
  )
  return(dataList)
}

#' Duplicate Datasets
#'
#' This function duplicates a data set a variable number of times. This function is to be used in
#' research situations.
#'
#' @param dataList A list of data set names that should be duplicated
#' @param duplicates The number of duplicates that are desired, defaults to 2
#' @return a data list of the duplicated datasets
#' @export
duplicate_datasets <- function(dataList, duplicates = 2)
{
  if(is.list(dataList)){
    names <- names(dataList)
    names <- rep(names, duplicates)
    dataList <- extract_datasets(names)
  }
  else
  {
    stop("A data list of datasets is required to duplicate datasets")
  }

  for (i in 1:length(names))
  {
    names[i] <- paste(i, names[i], sep="_")
    cols <- colnames(dataList[[i]])
    for (j in 1:length(cols))
    {
      cols[j] <- paste(i, cols[j], sep="_")
    }
    colnames(dataList[[i]]) <- cols
    names(dataList) <- names
  }
  return(dataList)
}



#' Permute the Rows of a Data List
#'
#' This function is designed for experimental analysis of algorithmic stability. It permutes the rows of a data list
#' so that further experimentation can be done with the data list to check the  algorithmic stability of a pipeline or
#' method
#'
#' @param dataList A list of data
#' @return The data list with the rows permuted
#' @export
permute_rows <- function(dataList)
{
  if(is.list(dataList)){
  for (i in (1:length(names(dataList)))) {
    dataList[[i]] <- dataList[[i]][sample(1:nrow(dataList[[i]])), ]
  }}
  else
  {
    stop("A data list of datasets is required to permute datasets rows")
  }
  return(dataList)
}


#' Permute the columns of a Data List
#'
#' This function is designed for experimental analysis of algorithmic stability. It permutes the columns of a data list
#' so that further experimentation can be done with the data list to check the  algorithmic stability of a pipeline or
#' method
#'
#' @param dataList A list of data
#' @return The data list with the columns permuted
#' @export
permute_columns <- function(dataList)
{
  if(is.list(dataList)){
  for (i in (1:length(names(dataList)))) {
    dataList[[i]] <- dataList[[i]][, sample(1:ncol(dataList[[i]]))]
  }}
  else
  {
    stop("A data list of datasets is required to permute datasets columns")
  }
  return(dataList)
}


#' Permute the Order fo the data sets
#'
#' This function is designed for experimental analysis of algorithmic stability.
#' This function permutes the order of the data set so that further experimentation can be done with the data list to check the  algorithmic stability of a pipeline or
#' method
#'
#' @param dataList a  list of data sets
#' @return A permuted data set list
#' @export
permute_dataset_order <- function(dataList)
{
  if(is.list(dataList)){
  for (i in (1:length(names(dataList)))) {
    permutedList <- dataList[sample(length(dataList))]
    return(permutedList)
  }}
  else
    {
      stop("A data list of datasets is required to permute datasets order")
    }
}


#' Modify Gene Counts
#'
#' This function takes the gene counts of a row of the takes the original gene counts of the datasets, multiplies them by 2 and then adds one to the count,
#' creating a transformation of 2x+1. This function is designed to test the algorithmic stability of a clustering algorithm
#'
#' @param dataList a  list of data sets
#' @return A data set list with modified gene counts
#' @export
modify_gene_counts <- function(dataList)
{
  if(is.list(dataList)){
    for (i in (1:length(names(dataList)))) {
      geneIndex <- sample.int(nrow(dataList[[i]]), 1)
      dataList[[i]][geneIndex,] <- 2 * dataList[[i]][geneIndex,] + 1
    }}
  else
  {
    stop("A list of datasets is required to modify gene counts.")
  }
  return(dataList)
}

#' Add duplicated Cells
#'
#' This function takes a random cell sample and adds its duplicate to the data set. This function is a test of algorithmic stability of a pipeline or method.
#'
#' @param dataList a  list of data sets
#' @return A data set list with a duplicated cell sample
#' @export
add_duplicate_cells <- function(dataList)
{
  if(is.list(dataList)){
    for (i in (1:length(names(dataList)))) {
      cellIndex <- sample.int(ncol(dataList[[i]]), 1)
      duplicatedCell <- dataList[[i]][, cellIndex]
      columnNames <- c(colnames(dataList[[i]]), paste(colnames(dataList[[i]])[cellIndex], '-dup', sep=''))
      dataList[[i]] <- cbind(dataList[[i]], duplicatedCell)
      colnames(dataList[[i]]) <- columnNames
    }}
  else
  {
    stop("A list of datasets is required to add duplicate cells.")
  }
  return(dataList)
}

#' Add Zero Variance Gene Counts
#'
#' This function adds a set of zero variance gene counts to the sample. This function is a test of algorithmic stability of a pipeline or a method.
#'
#' @param dataList a  list of data sets
#' @return A data set list with added zero variance gene counts
#' @export
add_zero_variance_gene_counts <- function(dataList)
{
  if(is.list(dataList)){
    for (i in (1:length(names(dataList)))) {
      zeroVarianceCounts <- rep(1, ncol(dataList[[i]]))
      dataList[[i]] <- rbind(dataList[[i]], zeroVarianceCounts)
    }}
  else
  {
    stop("A list of datasets is required to add zero variance gene counts.")
  }
  return(dataList)
}

#' Flips the Gene Counts
#'
#' This function multiplies the gene count values by -1, "flipping" them over an access. This function is a test of algorithmic stability of a pipeline or method.
#'
#' @param dataList a  list of data sets
#' @return A data set list with flipped gene counts
#' @export
flip_gene_counts <- function(dataList)
{
  if(is.list(dataList)){
    for (i in (1:length(names(dataList)))) {
      geneIndex <- sample.int(nrow(dataList[[i]]), 1)
      dataList[[i]][geneIndex,] <- (-1) * dataList[[i]][geneIndex,]
    }}
  else
  {
    stop("A list of datasets is required to filp the gene counts.")
  }
  return(dataList)
}


#' Run metamorphic tests
#'
#' This function runs the metamorphic test. It applies a perturbation and then clusters based on the method, transformation, and number of
#' clusters supplied.
#'
#' @param dataList A list of data sets
#' @param method The clustering method to be used, defaults to Kmeans
#' @param transformationType The type of data transformation to apply, defaults to log
#' @param seed the seed to be set for reproducibility, defaults to 1
#' @param numberClusters The number of cluster centers to use, defaults to 10
#' @return A metamorphic results table
#' @export
run_metamorphic_test <- function(dataList, perturbationType = 1, method = 'kmeans', transformationType = 'log', seed = 1, numberClusters = 10)
{
  if(is.list(dataList))
  {
    dataList = perturb_datasets(dataList, perturbationType)
    dataList = run_clustering_workflow(dataList, method, transformationType, seed, numberClusters)
    metamorphicReport = create_internal_cluster_validation_report(dataList, method)
  }
  else
  {
    stop("A list of datasets is required to run a metamorphic test")
  }
  return(metamorphicReport)
}

