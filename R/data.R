
#' Loads the stored data
#'
#' This function loads the data that is stored in the package for use in examples
#'
#' @export
data_load <- function(...)
{
  datasets <- readRDS(system.file("extdata", "gene_counts_v2.RDS", package = "benchmarking"))
  return(datasets)

}

#' Merge the Sparse Matricies
#'
#' This function
merge_sparse <- function(...)
{

  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()

  for (M in list(...)) {

    cnold <- colnames(M)
    rnold <- rownames(M)

    cnnew <- union(cnnew,cnold)
    rnnew <- union(rnnew,rnold)

    cindnew <- match(cnold,cnnew)
    rindnew <- match(rnold,rnnew)

    ind <- summary(M)
    i <- c(i,rindnew[ind[,1]])
    j <- c(j,cindnew[ind[,2]])
    x <- c(x,ind[,3])
  }
  Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
}


merge_datasets <- function(data.list, intersect=TRUE)
{
  data <- data.list[[1]]
  for (i in 2:length(data.list))
  {
    if(intersect) { data <- merge.sparse(data, data.list[[i]]) }
    else { data <- merge.Matrix(data, data.list[[i]], by.x = rownames(data), by.y=rownames(data.list[[i]]), all.x = TRUE, all.y = TRUE) }
  }
  return(list(data))
}

#' View Data
#'
#' This function displays a short summary of each data set in the data list provided. Can be used to gain quick insights into the
#' information that is stored in the data
#'
#' @param data.list A list of data sets that you would like to view
#' @export
view_data <- function(data.lsit)
{
  for (i in 1:length(data.list))
  {
    str(data.list[i])
  }
}

#' Extract Datasets
#'
#' This function extracts a set of Data sets for further use in the workflow. Each data set that is extracted will be a dgCMatrix
#' with rows representing genes and cells representing columns
#'
#' @param names A list of names representing the data sets you would like to extract
#' @return A data list of the matrices of the data sets that were extracted
#' @export
extract_datasets <- function(names)
{
  data.list <- gene_counts[names]
  for (name in names)
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


