
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

#' Merge the Sparse data sets
#'
#' This function merges the contents of two lists into a sparse matrix. It assumes that the data sets
#' are of the same format, with genes stored in rows and cells stored in columns. The lists are merged base on the
#' union of the row names and the union of the column names
#'
#' @param data.list A list of data that  you would like to merge
#' @return A list of the merged data sets
#' @export
merge_datasets <- function(data.list)
{
  cnnew <- character()
  rnnew <- character()
  x <- vector()
  i <- numeric()
  j <- numeric()

  for (M in data.list) {

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
#' TODO: FIX
#' @param data.list A list of data sets that you would like to view
#' @export
view_data <- function(idx)
{
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


#' Annotate Datasets
#'
#' This function annotates a created Seurat object. It edits the chunks and meta.data in
#' order to annotate the object for use in down stream analysis
#'
#' @param data A Seurat object that you want annotated
#' @return A annotated Seurat Object
#' @export
annotate_datasets <- function(data.list)
{
  data.list <- lapply(data.list, function(x)
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
    x@meta.data$SID <- as.numeric(as.factor(x@meta.data$ID))
    x@meta.data$orig.ident <- x@meta.data$CELL_LINE
    x <- x
  })
  return(data.list)
}
