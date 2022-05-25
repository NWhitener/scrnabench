#' Extract the Common Genes
#'
#' This function takes a data set and returns the genes that are "common" in the data set. It assumes that the row names
#' are the gene names.
#'
#'  @param data A scRNA-seq data object( Seurat) with the rows as genes and the cells as objects
#'  @return A list of the common genes in the data object
#'  @export
extract_common_genes <- function(data)
{
   common_gene_names <- Reduce(intersect, lapply(data, row.names))
   data.list <- lapply(data, function(x)
                { x[row.names(x) %in% common_gene_names,] })
   return(data.list)
}

#' Runs the GFICF Transformations
#'
#' This function completes the Gene Frequency Inverse Cell Frequency data transformation on a list of data. This step can be done to
#' create and use GFICF data for further downstream data analysis and computations. For more information on the GFICF calculation
#' please see https://github.com/dibbelab/gficf
#'
#' @param data.list A data list of Data sets in the Gene in Row and Cell in Columns format
#' @return A data list with GFICF scores of Each Gene Cell Relationship
#' @export
run_gficf <-function(data.list)
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- SeuratObject::CreateSeuratObject(counts = x)
    x$RNA@data <- gficf::gficf(M=x$RNA@counts, normalize=TRUE, storeRaw=FALSE)$gficf
    x$RNA@counts <- x$RNA@data[rownames(x$RNA@data),]
  })
  return(data.list)
}

#' Runs the Log Transformation
#'
#' This function completes the Log data transformation on a list of data. This step can be done to
#' create and use Log data for further downstream data analysis and computations. For more information on the Log calculation
#' please see https://cran.r-project.org/web/packages/Seurat/Seurat.pdf
#'
#' @param data.list A data list of Data sets in the Gene in Row and Cell in Columns format
#' @return A data list with Log scores of Each Gene Cell Relationship
#' @export
run_log <- function(data.list)
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- Seurat::NormalizeData(x)
  })
  return(data.list)
}

#' Preprocesses the a Data List
#'
#' This function "preprocesses" a data list by restricting a input list to a range between minimum and maximum of the Total MRNA counts
#' Total Gene Counts that are present in the list. This creates a filtered set of genes to use in experiments, removing extreme
#' outliers. These bounds are similar to the ?monocle? package. CITE PAPER
#'
#' @param data.list A list of data, with Genes in Rows and Cells in Columns
#' @return A data list of data that falls within the specified data limitations which acts as preprocesses
preprocess_data <- function(data.list)
{
  data.list <- lapply(X = data.list, function(x)
  {
    x <- Seurat::CreateSeuratObject(counts = x, min.cells = 3,min.features = 200)
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

  return(data.list)
}


# Make A test Dataset Example to make sure that we have a working preprocessing file/data file


