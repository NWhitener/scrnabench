
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




run_pca <- function(idx)
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- Seurat::RunPCA(x, verbose = FALSE)
  })
  return(data.list)
}
