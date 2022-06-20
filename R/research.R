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
