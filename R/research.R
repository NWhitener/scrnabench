
#Data Dependent
duplicate_datasets <- function(names, ndup)
{
  data.list <- list()
  if(length(names) == 2)
  {
    for (i in c(1: ndup))
    {
      x <- paste(names[1], "_dup", i, sep="")
      cells <- gene_counts[names[1]]
      data.list[x] <- cells
      x = paste(names[2], "_dup", i, sep="")
      cells <- gene_counts[names[2]]
      data.list[x] <- cells
    }
  }
  else
  {
    # there is a bug in this method
    for (i in c(1: ndup))
    {
      x <- paste(names[1], "_dup", i, sep="")
      cells <- gene_counts[names[1]]
      data.list[x] <- cells
    }
  }

  for (name in names(data.list))
  {
    data.list[[name]] <- as.matrix(data.list[[name]])
    data.list[[name]] <- as(data.list[[name]], "dgCMatrix")
    ix = Matrix::rowSums(data.list[[name]] != 0)
    data.list[[name]] = data.list[[name]][ix > 0,]
    ncols = length(colnames(data.list[[name]]))
    colnames(data.list[[name]]) = paste(name, seq(1:ncols), sep="-")
  }
  return(data.list)

}


#' Permute the Rows of a Data List
#'
#' This function is designed for experimental analysis of algorithmic stability. It permutes the rows of a data list
#' so that further experimentation can be done wit the data list to check the  algorithmic stability of a pipeline or
#' method
#'
#' @param data.list A data list that you would like to permute the rows of
#' @return The data list with the rows permuted
#' @export
permute_rows <- function(data.list)
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- x[sample(1:nrow(x)),]
  })
  return(data.list)
}


#' Permute the columns of a Data List
#'
#' This function is designed for experimental analysis of algorithmic stability. It permutes the columns of a data list
#' so that further experimentation can be done wit the data list to check the  algorithmic stability of a pipeline or
#' method
#'
#' @param data.list A data list that you would like to permute the columns of
#' @return The data list with the columns permuted
#' @export
permute_columns <- function(data.list)
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- x[,sample(1:ncol(x))]
  })
  return(data.list)
}



#' Run Principle Component Analysis
#'
#' This functions runs the Seurat RunPCA function on a list of data sets. This functions assumes that genes are in rows and
#' cells are in columns. The default parameters are used, and verbose is set to false
#'
#' @param data.list A data list of data sets
#' @return A data list with PCA completed on the features
#' @export
run_pca <- function(data.list)
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- Seurat::RunPCA(x, verbose = FALSE)
  })
  return(data.list)
}


#' Run Uniform Manifold Aproximation and Projection
#'
#' This functions runs the Seurat RunUMAP function on a list of data sets. This functions assumes that genes are in rows and
#' cells are in columns. The default reduction is set to "pca" with the first 30 dimensions being accepted
#'
#' @param data.list A data list of data sets
#' @return A data list with UMAP completed on the features
#' @export
run_umap <- function(data.list, reduction_choosen = "pca")
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- Seurat::RunUMAP(x,reduction = reduction_choosen, dims = 1:30)
  })
  return(data.list)
}

#' Complete Clustering Steps
#'
#' This functions runs the Seurat FindNeighbors and FindCLusters function on a list of data sets.
#' This functions assumes that genes are in rows and
#'  cells are in columns. The FindNeighbors reduction is set to "pca" by default and uses the first 30 dimensions
#'  The FindClusters resolution is set to 0.5
#'
#' @param data.list A data list of data sets
#' @return A data list with clustering completed completed on the features
#' @export
run_cluster <- function(data.list, reduction_choosen = "pca")
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- Seurat::FindNeighbors(x, reduction = reduction_choosen, dims = 1:30)
    x <- Seurat::FindClusters(x, resolution = 0.5)
  })
  return(data.list)
}
