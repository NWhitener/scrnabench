
#Data Dependent

#' Duplicate Datasets
#'
#' This function duplicates a dataset a variable number of times. This function is to be used in
#' research situations.
#'
#' @param names A list of dataset names that should be duplicated
#' @param ndups The number of duplicates that are desired
#' @return a dataList of the duplicated datasets
duplicate_datasets <- function(dataList, duplicates=2)
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
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- x[sample(1:nrow(x)),]
  })
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
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- x[,sample(1:ncol(x))]
  })
  return(dataList)
}



#' Run Principle Component Analysis
#'
#' This functions runs the Seurat RunPCA function on a list of data sets. This functions assumes that genes are in rows and
#' cells are in columns. The default parameters are used, and verbose is set to false
#'
#' @param dataList A data list of data sets
#' @return A data list with PCA completed on the features
#' @export
run_pca <- function(dataList)
{
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- Seurat::RunPCA(x, verbose = FALSE)
  })
  return(dataList)
}


#' Run Uniform Manifold Aproximation and Projection
#'
#' This functions runs the Seurat RunUMAP function on a list of data sets. This functions assumes that genes are in rows and
#' cells are in columns. The default reduction is set to "pca" with the first 30 dimensions being accepted
#'
#' @param dataList A data list of data sets
#' @return A data list with UMAP completed on the features
#' @export
run_umap <- function(dataList, reduction_choosen = "pca")
{
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- Seurat::RunUMAP(x,reduction = reduction_choosen, dims = 1:30)
  })
  return(dataList)
}

#' Complete Clustering Steps
#'
#' This functions runs the Seurat FindNeighbors and FindCLusters function on a list of data sets.
#' This functions assumes that genes are in rows and
#'  cells are in columns. The FindNeighbors reduction is set to "pca" by default and uses the first 30 dimensions
#'  The FindClusters resolution is set to 0.5
#'
#' @param dataList A data list of data sets
#' @return A data list with clustering completed completed on the features
#' @export
run_cluster <- function(dataList, reduction_choosen = "pca")
{
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- Seurat::FindNeighbors(x, reduction = reduction_choosen, dims = 1:30)
    x <- Seurat::FindClusters(x, resolution = 0.5)
  })
  return(dataList)
}



permute_dataset_order <- function(dataList)
{
  permuted.list <- dataList[sample(length(dataList))]
  return(permuted.list)
}
