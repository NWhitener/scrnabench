#' Extract the Common Genes
#'
#' This function takes a data set and returns the genes that are "common" in the data set. It assumes that the row names
#' are the gene names.
#'
#' @param data A scRNA-seq data object( Seurat) with the rows as genes and the cells as objects
#' @return A list of the common genes in the data object
#' @export
extract_common_genes <- function(dataList)
{
  ###Error Handling/warning
  if(is.list(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to extract common genes")
  }
  if(length(dataList) == 1)
  {
    warning("Only one dataset provided, returning the original dataset")
    return(dataList)
  }
  common_gene_names <- Reduce(intersect, lapply(dataList, row.names))
  dataList <- lapply(dataList, function(x)
  { x[row.names(x) %in% common_gene_names,] })
  return(dataList)
}

#' Runs the GFICF Transformations
#'
#' This function completes the Gene Frequency Inverse Cell Frequency data transformation on a list of data. This step can be done to
#' create and use GFICF data for further downstream data analysis and computations. For more information on the GFICF calculation
#' please see https://github.com/dibbelab/gficf
#'
#' @param dataList A data list of Data sets in the Gene in Row and Cell in Columns format
#' @return A data list with GFICF scores of Each Gene Cell Relationship
#' @export
run_gficf <-function(dataList)
{

  if(is.list(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to apply GFICF to datasets")
  }
    for (i in range(1:length(names(data.list))))
    {
      data.list[[i]] <- Seurat::CreateSeuratObject(counts = data.list[[i]])
      data.list[[i]]$RNA@data <- gficf::gficf(M=data.list[[i]]$RNA@counts,
                                              normalize=TRUE, storeRaw=FALSE)$gficf
      data.list[[i]]$RNA@counts <- data.list[[i]]$RNA@data[rownames(data.list[[i]]$RNA@data),]
    }
    return(data.list)
}

#' Runs the Log Transformation
#'
#' This function completes the Log data transformation on a list of data. This step can be done to
#' create and use Log data for further downstream data analysis and computations. For more information on the Log calculation
#' please see https://cran.r-project.org/web/packages/Seurat/Seurat.pdf
#'
#' @param dataList A data list of Data sets in the Gene in Row and Cell in Columns format
#' @return A data list with Log scores of Each Gene Cell Relationship
#' @export
run_log <- function(dataList)
{
  if(is.list(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to apply the Log transformations to datasets")
  }
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- Seurat::NormalizeData(x)
  })
  return(dataList)
}

#' Preprocesses the a Data List
#'
#' This function "preprocesses" a data list by restricting a input list to a range between minimum and maximum of the Total MRNA counts
#' Total Gene Counts that are present in the list. This creates a filtered set of genes to use in experiments, removing extreme
#' outliers. These bounds are similar to the ?monocle? package. CITE PAPER
#'
#' @param dataList A list of data, with Genes in Rows and Cells in Columns
#' @return A data list of data that falls within the specified data limitations which acts as preprocesses
#' @export
preprocess <- function(dataList)
{
  if(is.list(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to preprocess datasets")
  }

  dataList <- lapply(X = dataList, function(x)
  {
    x <- Seurat::CreateSeuratObject(counts = x, min.cells = 3, min.features = 200)
    x[["percent.mt"]] <- Seurat::PercentageFeatureSet(x, pattern = "^MT-")
    Total_mRNAs <- x[["nCount_RNA"]]$nCount_RNA
    mupper_bound <- 10^(mean(log10(Total_mRNAs)) + 2*sd(log10(Total_mRNAs)))
    mlower_bound <- 10^(mean(log10(Total_mRNAs)) - 2*sd(log10(Total_mRNAs)))
    Total_Genes <- x[["nFeature_RNA"]]$nFeature_RNA
    gupper_bound <- 10^(mean(log10(Total_Genes)) + 2*sd(log10(Total_Genes)))
    glower_bound <- 10^(mean(log10(Total_Genes)) - 2*sd(log10(Total_Genes)))
    x <- subset(x = x, subset = nFeature_RNA > glower_bound & nFeature_RNA < gupper_bound &
                  nCount_RNA > mlower_bound & nCount_RNA < mupper_bound & percent.mt < 10)
  })

  return(dataList)
}

#' Select Highly Variable Genes
#'
#' This function uses the FindVariablesFeatures of Seurat to find the Highly variable genes of a data set. This function
#' assumes that genes are in rows and cells are in columns. The "vst" selection method is used with 2000 nfeatures.
#'
#' @param dataList A list of data sets that you wouuld like the variable genes of
#' @return A data list of that includes the variable genes of the original data list
#' @export
select_hvg <- function(dataList)
{
  if(is.list(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to select the Highly Variable Genes of a dataset")
  }
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  return(dataList)
}

#' Scales the Data
#'
#' This function uses the ScaleData of Seurat to scale the data set. This function
#' assumes that genes are in rows and cells are in columns.
#'
#' @param dataList A list of data sets that you would like to scale
#' @return A data list of that includes the scaled data
#' @export
scale_data <- function(dataList)
{

  if(is.list(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to scale datasets")
  }
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- Seurat::ScaleData(x)
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

#' Run tSNE
#'
#' This functions runs the Seurat RunTSNE function on a list of data sets. This functions assumes that genes are in rows and
#' cells are in columns. The default reduction is set to "pca" with the first 30 dimensions being accepted
#'
#' @param dataList A data list of data sets
#' @return A data list with tSNE completed on the features
#' @export
run_tsne <- function(dataList, reduction_choosen = 'pca')
{
  if(is.list(dataList)== "FALSE")
  {
    stop("A data list of datasets is required to apply the tSNE reduction to datasets")
  }
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- Seurat::RunTSNE(x, reduction = reduction_choosen, dims = 1:30)
  })
  return(dataList)
}
