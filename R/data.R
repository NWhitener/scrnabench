
#' Loads the stored data
#'
#' This function loads the data that is stored in the package for use in examples
#'
#' @export
data_load <- function()
{

  gene_counts <<- readRDS(system.file("extdata", "gene_counts_v6.RDS", package = "benchmarking", mustWork = TRUE))
  datasets <<- ls(gene_counts)
  return(datasets)

}

#' Merge the Sparse Matrices
#'
#' This function merges the contents of two lists into a sparse matrix. It assumes that the data sets
#' are of the same format, with genes stored in rows and cells stored in columns. The lists are merged base on the
#' union of the row names and the union of the column names
#'
#' TODO: FIX arguement passing
#'
#' @export
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

#' Merges Data sets
#'
#' This function merges a series of data lists into one data matrix, based on the intersection of the individuals. If there is intersection
#' the matrix is merged as sparse to handle zeros, otherwise it is merged as a matrix. This function assumes that the data is
#' in the format genes in rows and cells in columns.
#'
#' @param data.list A list of data set to merge
#' @param intersect A Boolean parameter, that decides if the output should be a sparse matrix or a matrix
#' @return A data list of the merged data
#' @export
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
#' TODO: FIX
#' @param data.list A list of data sets that you would like to view
#' @export
view_data <- function(data)
{
  for (i in 1:length(data))
  {
    str(extract_datasets(data[i]))

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


