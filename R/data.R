

#' Downloads the data
#'
#' This function downloads the data that is used in the package for use in examples, and automatically loads the data for
#' usage. Run this once when the data is not download, if the data is downloaded use the load_data function
#'
#' @import inborutils
#' @param path The path of the download location to use
#' @return A object with the names of the data sets that were downloaded and loaded
#' @export
download_data <- function(path = '.')
{

  inborutils::download_zenodo(doi = "10.5281/zenodo.6617997", path = path)
  geneCounts <<- readRDS(file = paste(path, "/gene_counts.RDS", sep = ''))
  datasets <- ls(geneCounts)

  return(datasets)
}


#' Loads the data
#'
#' This function loads the data that is used in the package for use in workflows, functions, and examples
#' @param demo Use the full dataset or the small demo data set, defaulted to use the full data sets
#' @param path The path of the load location to use, if the full datasets is to be loaded
#' @return A object with the names of the data sets that were loaded
#' @export
load_data <- function(demo = FALSE, path = '.')
  {

    path_full = paste(path, "/gene_counts.RDS", sep = '')
    if(demo)
    {
      geneCounts  <<- readRDS(system.file("extdata", "gene_counts_demo.RDS", package = "scrnabench", mustWork = TRUE))
    }
    else if(file.exists(path_full))
    {
      geneCounts <<- readRDS(file = path_full)
    }
    else
    {
      stop("Please Provide a Valid File Path, or use the download_data() function to download the data")
    }
    datasets <- ls(geneCounts)

    return(datasets)
  }

#' Merge the Sparse data sets
#'
#' This function merges the contents of two lists into a sparse matrix. It assumes that the data sets
#' are of the same format, with genes stored in rows and cells stored in columns. The lists are merged base on the
#' union of the row names and the union of the column names
#'
#' @import Matrix
#' @param dataList A list of data that  you would like to merge
#' @return A list of the merged data sets
#' @export
merge_datasets <- function(dataList)
{

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

  for (data in dataList)
  {

    cnold <- colnames(data)
    rnold <- rownames(data)

    cnnew <- union(cnnew,cnold)
    rnnew <- union(rnnew,rnold)

    cindnew <- match(cnold,cnnew)
    rindnew <- match(rnold,rnnew)

    ind <- Matrix::summary(data)
    i <- c(i,rindnew[ind[,1]])
    j <- c(j,cindnew[ind[,2]])
    x <- c(x,ind[,3])
  }
  merged <- list(Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(length(rnnew),length(cnnew)),dimnames=list(rnnew,cnnew)))
  names(merged) <- c("Integrated")

  return(merged)
}

#' Describe Datasets
#'
#' This function displays a short summary of each data set in the data list provided. Can be used to gain quick insights into the
#' information that is stored in the data
#'
#' @param dataList A list of names of data sets that you would like to view
#' @export
describe_datasets <- function(dataList)
{

  # Show the dataset name, number of genes, number of cells

  summaryTable = NULL
  for (i in 1:length(names(dataList)))
  {
    datasetNames = names(dataList[i])
    numberGenes <- as.numeric(nrow(dataList[[i]]))
    numberCells <- as.numeric(ncol(dataList[[i]]))
    summaryTable <- rbind(summaryTable, cbind(datasetNames, numberGenes, numberCells))
  }
  colnames(summaryTable) <- c('Dataset Name', 'Number of Genes', 'Number of Cells')
  print(as.data.frame(summaryTable))
}



#' Extract Datasets
#'
#' This function extracts a set of Data sets for further use in the workflows. Each data set that is extracted will be a dgCMatrix
#' with rows representing genes and cells representing columns
#'
#' @import Matrix
#' @param names A list of names representing the data sets you would like to extract
#' @return A data list of the matrices of the data sets that were extracted
#' @export
extract_datasets <- function(names)
{
  ##Check to see geneCounts exists

if(exists('geneCounts')){
  if(!is.vector(names))
  {
    stop("A list of names of the datasets is required to extract datasets")
  }
  dataList <- geneCounts[names]
  for (name in names)
  {
    dataList[[name]] <- as.matrix(dataList[[name]])
    dataList[[name]] <- methods::as(dataList[[name]], "dgCMatrix")
    ix = Matrix::rowSums(dataList[[name]] != 0)
    dataList[[name]] = dataList[[name]][ix > 0,]
    ncols = length(colnames(dataList[[name]]))
    colnames(dataList[[name]]) = paste(name, seq(1:ncols), sep="-")
  }
  return(dataList)
}
  else{
    stop("Please run load_data() or download_data() first")
  }
}


#' Annotate Data sets
#'
#' This function annotates a created Seurat object. It edits the chunks and meta.data in
#' order to annotate the object for use in down stream analysis
#'
#' @param dataList A list of Seurat object that you want annotated
#' @return A annotated Seurat Object
#' @export
annotate_datasets <- function(dataList)
{
  if(is.list(dataList))
  {
    dataList <- lapply(dataList, function(x)
    {
      cols = colnames(x)
      for (i in (1:length(cols)))
      {
        col <- cols[i]
        col <- stringr::str_split(col, "-")[[1]][1]
        x@meta.data$ID[i] <- col

        chunks <- stringr::str_split(col, "_")
        x@meta.data$TECHNOLOGY[i] <- chunks[[1]][1]

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
  }

  else{
    stop("A data list of datasets is required to annotate datasets")
  }
  return(dataList)
}


#' Add Clustering Annotations
#'
#' This function will add clustering annotations (i.e memberships) to the metadata of the Seurat Object of the matching data.
#' Any clustering method can be used to obtain these memberships. The annotation name will include the clustering type
#' and the cells cluster ID.
#'
#' @param dataList A list of data sets to add the clustering annotation too
#' @param memberships A list of the cluster memberships of each cell, used to add the annotation
#' @param clusteringName The name of the clustering method used to obtain the clustering membership
#' @return A list of the data sets with the clustering annotations added
#' @export
add_clustering_annotation <- function(dataList, memberships, clusteringName)
{
  if(is.list(dataList) & is.list(memberships))
  {
    for (i in (1: length(names(dataList))))
    {
      dataList[[i]][[clusteringName]] <- memberships[i]
    }
  }
  else
  {
    stop("A list of datasets and memberships is required to annotate datasets")
  }
  return(dataList)
}
