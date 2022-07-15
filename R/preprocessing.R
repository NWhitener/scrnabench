#' Extract the Common Genes
#'
#' This function takes a data set and returns the genes that are "common" in the data set. It assumes that the row names
#' are the gene names.
#'
#' @param dataList A List of scRNA-seq data object (Seurat) with the rows as genes and the cells as objects
#' @return A list of the common genes in the data object
#' @export
extract_common_genes <- function(dataList)
{

  if(is.list(dataList))
  {
    if(length(dataList) == 1)
    {
      warning("Only one dataset provided, returning the original dataset")
    }
    else
    {
      common_gene_names <- Reduce(intersect, lapply(dataList, row.names))
      dataList <- lapply(dataList, function(x)
      {
        x[row.names(x) %in% common_gene_names,] })
    }
  }
  else
  {
    stop("A data list of datasets is required to extract common genes")
  }

  return(dataList)
}

#' Runs the tfidf Transformations
#'
#' This function completes the Term Frequency Inverse Document Frequency data transformation on a list of data. This is step is complete by
#' using the text2vec protocol for completing TF-IDF on the cell data, creating TF-IDF Scores
#'
#' @param dataList A data list of Data sets in the Gene in Row and Cell in Columns format
#' @return A data list with GFICF scores of Each Gene Cell Relationship
#' @export
run_tfidf<-function(dataList)
{
  if(is.list(dataList))
  {
    model <- text2vec::TfIdf$new(smooth_idf = TRUE, norm = c('l1'), sublinear_tf = TRUE)
    for (i in (1:length(names(dataList))))
    {
      tfidf <- model$fit_transform(t(dataList[[i]]$RNA@data))
      dataList[[i]]$RNA@data <- t(tfidf)
      dataList[[i]] <- Seurat::SetAssayData(object=dataList[[i]], slot = 'scale.data', new.data = as.matrix(dataList[[i]]$RNA@data))
    }
  }
  else
  {
    stop("A data list of datasets is required to apply TFIDF to datasets")
  }

  return(dataList)
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
  if(is.list(dataList))
  {
    for (i in (1:length(names(dataList))))
    {
      dataList[[i]]<- Seurat::NormalizeData(dataList[[i]])
    }
  }
  else
  {
    stop("A data list of datasets is required to apply the Log transformations to datasets")
  }

  return(dataList)
}

#' Filter the a Data List
#'
#' This function "filters" a data list by restricting a input list to a range between minimum and maximum of the Total MRNA counts
#' Total Gene Counts that are present in the list. This creates a filtered set of genes to use in experiments, removing extreme
#' outliers and making sure that all cells and genes are potentially informative.
#'
#' @param dataList A list of data, with Genes in Rows and Cells in Columns
#' @return A data list of data that falls within the specified data limitations which acts as preprocesses
#' @export
filter_data <- function(dataList)
{

  nCount_RNA = NULL
  nFeature_RNA = NULL
  percent.mt = NULL
  if(is.list(dataList))
  {
    for (i in (1:length(names(dataList))))
    {
      dataList[[i]]<- Seurat::CreateSeuratObject(counts = dataList[[i]], min.cells = 3, min.features = 200)
      dataList[[i]][["percent.mt"]] <- Seurat::PercentageFeatureSet(dataList[[i]], pattern = "^MT-")
      Total_mRNAs <- dataList[[i]][["nCount_RNA"]]$nCount_RNA
      mupper_bound <- 10^(mean(log10(Total_mRNAs)) + 2*stats::sd(log10(Total_mRNAs)))
      mlower_bound <- 10^(mean(log10(Total_mRNAs)) - 2*stats::sd(log10(Total_mRNAs)))
      Total_Genes <- dataList[[i]][["nFeature_RNA"]]$nFeature_RNA
      gupper_bound <- 10^(mean(log10(Total_Genes)) + 2*stats::sd(log10(Total_Genes)))
      glower_bound <- 10^(mean(log10(Total_Genes)) - 2*stats::sd(log10(Total_Genes)))
      dataList[[i]] <- subset(x = dataList[[i]], subset = nFeature_RNA > glower_bound & nFeature_RNA < gupper_bound &
                                nCount_RNA > mlower_bound & nCount_RNA < mupper_bound & percent.mt < 10)
    }
  }
  else{
    stop("A data list of datasets is required to preprocess datasets")
  }

  return(dataList)
}

#' Select Highly Variable Genes
#'
#' This function uses the FindVariablesFeatures of Seurat to find the Highly variable genes of a data set. This function
#' assumes that genes are in rows and cells are in columns. The "vst" selection method is used with 2000 genes as a default.
#'
#' @param dataList A list of data sets for finding HVG's
#' @param numGenes The number fo HVG's to be found, defaults to 2000
#' @return A data list of that includes the variable genes of the original data list
#' @export
select_hvg <- function(dataList, numGenes = 2000)
{
  if(is.list(dataList))
  {

    for (i in (1:length(names(dataList))))
    {
      dataList[[i]] <- Seurat::FindVariableFeatures(dataList[[i]], selection.method = "vst", nfeatures = numGenes)
    }
  }
  else
  {
    stop("A data list of datasets is required to select the Highly Variable Genes of a dataset")
  }

  return(dataList)
}

#' Scales the Data
#'
#' This function uses the ScaleData of Seurat to scale the data set. Scaling centers the data to a mean of zero,
#' and scales to a standard deviation of 1. This function assumes that genes are in rows and cells are in columns.
#'
#' @param dataList A list of data sets that you would like to scale
#' @return A data list of that includes the scaled data
#' @export
scale_data <- function(dataList)
{

  if(is.list(dataList))
  {
    for (i in (1:length(names(dataList))))
    {
      dataList[[i]] <- Seurat::ScaleData(dataList[[i]])
    }
  }
  else
  {
    stop("A data list of datasets is required to scale datasets")
  }
  return(dataList)
}


#' Run Principle Component Analysis
#'
#' This functions runs the Seurat RunPCA function on a list of data sets. This functions assumes that genes are in rows and
#' cells are in columns. Verbose is set to false, the Features are set to be the Variable Features of the data set, the
#' default number of components is 10, and the PCA is run on only scaled data
#'
#' @param dataList A list of data sets to be reduced by PCA
#' @param numComponents The number fo components to be found, with a default of 10
#' @return A data list with PCA completed on the features
#' @export
run_pca <- function(dataList, numComponents = 10)
{
  if(is.list(dataList))
  {
    for (i in (1:length(names(dataList))))
    {
      dataList[[i]] <- Seurat::RunPCA(dataList[[i]], features = Seurat::VariableFeatures(object = dataList[[i]]),
                                      npcs = numComponents, verbose = FALSE)
      dataList[[i]] <- Seurat::DietSeurat(dataList[[i]], data = T,
                                          counts = F, scale.data = F, features = Seurat::VariableFeatures(object = dataList[[i]]),
                                          dimreducs = names(dataList[[i]]@reductions))
    }
  }
  else
  {
    stop("A data list of datasets is required to scale datasets")
  }
  return(dataList)
}



#' Run UMAP
#'
#' This functions runs the Seurat RunUMAP function on a list of data sets. This functions assumes that genes are in rows and
#' cells are in columns. The default reduction is set to "pca" with the first 10 dimensions being accepted
#'
#' @param dataList A list of data sets to be reduced by UMAP
#' @param reductionType The type of reduction that the UMAP calculation should be based on, defaulted to PCA
#' @param numDimensions The number of dimensions to be used in the UMAP calculation, default is 10
#' @return A data list with UMAP completed on the features
#' @export
run_umap <- function(dataList, reductionType = 'pca', numDimensions = 10)
{

  if(is.list(dataList))
  {
    for (i in (1:length(names(dataList))))
    {
      dataList[[i]] <- Seurat::RunUMAP(dataList[[i]], reduction = reductionType, dims = 1:numDimensions)
    }
  }
  else
  {
    stop("A data list of datasets is required to reduce datasets using UMAP")
  }

  return(dataList)
}

#' Run tSNE
#'
#' This functions runs the Seurat RunTSNE function on a list of data sets. This functions assumes that genes are in rows and
#' cells are in columns. The default reduction is set to "pca" with the first 10 dimensions being accepted. The perplexity
#' is set to the square root of the number of cells, to accomodate for small datasets, and check_duplicated is set
#' to be false
#'
#' @param dataList A list of data sets to be reduced by tSNE
#' @param reductionType The type of reduction that the tSNE calculation should be based on, defaulted to PCA
#' @param numDimensions The number of dimensions to be used in the tSNE calculation, default is 10
#' @return A data list with tSNE completed on the features
#' @export
run_tsne <- function(dataList, reductionType = 'pca', numDimensions = 10)
{
  if(is.list(dataList))
  {
    for (i in (1:length(names(dataList))))
    {
      perplexityValue <- sqrt(ncol(dataList[[i]]))
      dataList[[i]] <- Seurat::RunTSNE(dataList[[i]], reduction = reductionType, dims = 1:numDimensions,
                                       perplexity = perplexityValue, check_duplicates = FALSE)
    }
  }
  else
  {
    stop("A data list of datasets is required to reduce datasets using TSNE")
  }
  return(dataList)
}
