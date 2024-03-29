% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate.R
\name{run_silhouette}
\alias{run_silhouette}
\title{Run Silhouette}
\usage{
run_silhouette(dataList, method = "kmeans", reductionType = "pca")
}
\arguments{
\item{dataList}{a data list of the data with clustering assignments}

\item{method}{The clustering method used, defaults to Kmeans}

\item{reductionType}{The dimensionality reduction type of the data to find the calculated scores,
defaults to PCA}
}
\value{
A table of the Silhouette Scores, with the data set name and the score
}
\description{
This function calculates the average silhouette width of the clustering assignments made on the data sets provided.
Uses the clusterCrit package implementation.
}
