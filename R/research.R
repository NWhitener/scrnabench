

#' Duplicate Datasets
#'
#' This function duplicates a dataset a variable number of times. This function is to be used in
#' research situations.
#'
#' @param names A list of dataset names that should be duplicated
#' @param ndups The number of duplicates that are desired
#' @return a dataList of the duplicated datasets
duplicate_datasets <- function(dataList, duplicates = 2)
{
  names <- names(dataList)
  names <- rep(names, duplicates)
  dataList <- extract_datasets(names)
  return(dataList)
}


#' Permute the Rows of a Data List
#'
#' This function is designed for experimental analysis of algorithmic stability. It permutes the rows of a data list
#' so that further experimentation can be done wit the data list to check the  algorithmic stability of a pipeline or
#' method
#'
#' @param dataList A data list that you would like to permute the rows of
#' @return The data list with the rows permuted
#' @export
permute_rows <- function(dataList)
{
  for (i in (1:length(names(dataList)))) {
    dataList[[i]] <- dataList[[i]][sample(1:nrow(dataList[[i]])), ]
  }
  return(dataList)
}


#' Permute the columns of a Data List
#'
#' This function is designed for experimental analysis of algorithmic stability. It permutes the columns of a data list
#' so that further experimentation can be done wit the data list to check the  algorithmic stability of a pipeline or
#' method
#'
#' @param dataList A data list that you would like to permute the columns of
#' @return The data list with the columns permuted
#' @export
permute_columns <- function(dataList)
{
  for (i in (1:length(names(dataList)))) {
    dataList[[i]] <- dataList[[i]][, sample(1:ncol(dataList[[i]]))]
  }
  return(dataList)
}


#' Permute the Order fo the data sets
#'
#' This function permutes the order of the data set
#'
#' @param dataList a data list of data sets to permute the order of
#' @return A permuted dataset list
#' @export
permute_dataset_order <- function(dataList)
{
  for (i in (1:length(names(dataList)))) {
    permutedList <- dataList[sample(length(dataList))]
    return(permutedList)
  }
}
