
#' Loads the stored data
#'
#' This function loads the data that is stored in the package for use in examples
#' @param full Download the full dataset or use the small demo dataset
#' @param path The path of the download location to use
#' @export
data_load <- function(full = FALSE, path = '.')
{
  if(full == TRUE)
  {
    inborutils::download_zenodo(doi = "10.5281/zenodo.6617997", path = path)
    gene_counts <<- readRDS(file = paste(path, "/gene_counts.RDS", sep = ''))
  }
  else{
  gene_counts <<- readRDS(system.file("extdata", "gene_counts_v6.RDS", package = "benchmarking", mustWork = TRUE))
  }
  datasets <<- ls(gene_counts)
  return(datasets)

}

#' Merge the Sparse data sets
#'
#' This function merges the contents of two lists into a sparse matrix. It assumes that the data sets
#' are of the same format, with genes stored in rows and cells stored in columns. The lists are merged base on the
#' union of the row names and the union of the column names
#'
#' @param dataList A list of data that  you would like to merge
#' @return A list of the merged data sets
#' @export
merge_datasets <- function(dataList)
{

  if(is.list(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to merge datasets")
  }
  if(length(dataList) == 1)
  {
    warning("Only one dataset provided, returning the original dataset")
    return(dataList)
  }
  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()

  for (M in dataList) {

    cnold <- colnames(M)
    rnold <- rownames(M)

    cnnew <- union(cnnew,cnold)
    rnnew <- union(rnnew,rnold)

    cindnew <- match(cnold,cnnew)
    rindnew <- match(rnold,rnnew)

    ind <- Matrix::summary(M)
    i <- c(i,rindnew[ind[,1]])
    j <- c(j,cindnew[ind[,2]])
    x <- c(x,ind[,3])
  }
  merged <- Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew))
  return(list(merged))
}




#' View Data
#'
#' This function displays a short summary of each data set in the data list provided. Can be used to gain quick insights into the
#' information that is stored in the data
#'
#' @param dataList A list of data sets that you would like to view
#' @export
view_data <- function(idx)
{

  if(!is.vector(idx))
  {
    stop("A list of names of the datasets is required to view datasets")
  }
  for (i in 1:length(idx))
  {
    str(extract_datasets(idx[i]))

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
  if(!is.vector(names))
  {
    stop("A list of names of the datasets is required to extract datasets")
  }
  dataList <- gene_counts[names]
  for (name in names)
  {
    dataList[[name]] <- as.matrix(dataList[[name]])
    dataList[[name]] <- as(dataList[[name]], "dgCMatrix")
    ix = Matrix::rowSums(dataList[[name]] != 0)
    dataList[[name]] = dataList[[name]][ix > 0,]
    ncols = length(colnames(dataList[[name]]))
    colnames(dataList[[name]]) = paste(name, seq(1:ncols), sep="-")
  }
  return(dataList)
}


#' Annotate Datasets
#'
#' This function annotates a created Seurat object. It edits the chunks and meta.data in
#' order to annotate the object for use in down stream analysis
#'
#' @param data A Seurat object that you want annotated
#' @return A annotated Seurat Object
#' @export
annotate_datasets <- function(dataList)
{
  if(is.list(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to annotate datasets")
  }
  dataList <- lapply(dataList, function(x)
  {
    cols = colnames(x)
    for (i in (1:length(cols)))
    {
      col <- cols[i]
      col <- stringr::str_split(col, "-")[[1]][1]
      x@meta.data$ID[i] <- col

      chunks <- stringr::str_split(col, "_")
      x@meta.data$TECHNLOGY[i] <- chunks[[1]][1]

      if(chunks[[1]][2] == "PE" | chunks[[1]][2] == "SE")
      {
        x@meta.data$CENTER[i] <- "TBU"
      }
      else
      {
        x@meta.data$CENTER[i] <- chunks[[1]][2]
      }

      if(chunks[[1]][3] == "M" | chunks[[1]][3] == "HT")
      {
        x@meta.data$CELL_LINE[i] <- ifelse(chunks[[1]][4] == "A", 'HCC1395', 'HCC1395BL')
        x@meta.data$PREPROCESS[i] <- chunks[[1]][5]
      }
      else
      {
        x@meta.data$CELL_LINE[i] <- ifelse(chunks[[1]][3] == "A", 'HCC1395', 'HCC1395BL')
        x@meta.data$PREPROCESS[i] <- chunks[[1]][4]
      }
    }
    x@meta.data$SID <- rep(1, length(x@meta.data$ID))
    x@meta.data$orig.ident <- x@meta.data$CELL_LINE
    x <- x
  })

  for (i in range(1, length(dataList)))
  {
    dataList[[i]]@meta.data$SID = rep(i, length(dataList[[i]]@meta.data$SID))
  }

  return(dataList)
}

