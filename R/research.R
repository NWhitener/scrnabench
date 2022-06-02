
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


#' Permute the Order fo the data sets
#'
#' This function permutes the order of the data set
#'
#' @param dataList a data list of data sets to permute the order of
#' @return A permuted dataset list
#' @export
permute_dataset_order <- function(dataList)
{
  dataList <- lapply(X = dataList, FUN = function(x) {
  permutedList <- dataList[sample(length(dataList))]
  return(permutedList)
  })
}


#' Run Kmeans
#'
#' This function runs Kmeans Clustering. It uses the simple Kmeans algorithm with
#' a max nunber of iterations set to be 100.
#'
#' @param dataList A list of datasets to apply the Kmeans algorithm too
#' @param k The number of cluster centers to use, defaults to 10
#' @return The dataList with kmeans clusters stored in the meta data
#' @export
run_kmeans <- function(dataList, k=10)
{
  if(is.vector(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to apply the cluster the datasets")
  }

  dataList <- lapply(X = dataList, FUN = function(x){
    clustData = x@reductions[["pca"]]@cell.embeddings
    x@meta.data$kmeans_cluster<- stats::kmeans(clustData,  k, iter.max = 100)$cluster
    x <- x
    print(x@meta.data)
  })

  return(dataList)
}
