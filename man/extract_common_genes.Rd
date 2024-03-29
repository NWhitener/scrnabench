% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{extract_common_genes}
\alias{extract_common_genes}
\title{Extract the Common Genes}
\usage{
extract_common_genes(dataList)
}
\arguments{
\item{dataList}{A List of scRNA-seq data object (Seurat) with the rows as genes and the cells as objects}
}
\value{
A list of the common genes in the data object
}
\description{
This function takes a data set and returns the genes that are "common" in the data set. It assumes that the row names
are the gene names.
}
