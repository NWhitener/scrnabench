#' Visualization method for the Clustering Methods
#'
#' Produces a plot of the clustering, based on the dimensionality reduction method chosen. The plot is a Seurat DimPlot
#'
#' @param data An individual data object that you want to visualize
#' @param reductionTypeThe dimensionality reduction method chosen, defaults to PCA
#' @param method The cluster method to use, defaults to Kmeans
#' @export
plot_clusters <- function(data, reductionType = 'pca', method = 'kmeans')
{
  annotationField <- toupper(paste(method, '_cluster_', reductionType, sep=''))
  return(Seurat::DimPlot(data, reduction = reductionType, group.by = annotationField))
}

