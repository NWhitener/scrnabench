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
  return(dataList)
}


#' Permute the Rows of a Data List
#'
#' This function is designed for experimental analysis of algorithmic stability. It permutes the rows of a data list
#' so that further experimentation can be done wit the data list to check the  algorithmic stability of a pipeline or
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
#' so that further experimentation can be done wit the data list to check the  algorithmic stability of a pipeline or
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
#' This function permutes the order of the data set
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
#' This function takes the gene counts of the takes the original gene counts of the datasets, multiplies them by 2 and then adds one to the count,
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
#' This function takes a random cell sample and adds its duplicate to the data set. This function is a test of algorithmic stability
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
#' This function adds a set of zero variance gene counts to the sample. This function is a test of algorithmic stability
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
#' This function multiplies the gene count values by -1, "flipping" them over an access. This function is a test of algorithmic stability
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
