% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/research.R
\name{modify_gene_counts}
\alias{modify_gene_counts}
\title{Modify Gene Counts}
\usage{
modify_gene_counts(dataList)
}
\arguments{
\item{dataList}{a  list of data sets}
}
\value{
A data set list with modified gene counts
}
\description{
This function takes the gene counts of a row of the takes the original gene counts of the datasets, multiplies them by 2 and then adds one to the count,
creating a transformation of 2x+1. This function is designed to test the algorithmic stability of a clustering algorithm
}
