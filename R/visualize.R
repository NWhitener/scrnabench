#' Visualization method for the Clustering Methods
#'
#' Produces a plot of the clustering, based on the dimensionality reduction method chosen. The plot is a Seurat DimPlot.
#'
#' @param data An individual data object that you want to visualize
#' @param reductionType The dimensionality reduction method chosen, defaults to PCA
#' @param method The cluster method to use, defaults to Kmeans
#' @param transformationType The transformation type you would like to see appear in the title
#' @export
plot_clusters <- function(data,method = 'kmeans', reductionType = 'pca', transformationType = 'log')
{
  annotationField <- toupper(paste(method, '_cluster_', reductionType, sep=''))
  titlePlot =  paste(toupper(method), "_", toupper(reductionType), "_", toupper(transformationType), sep = "")
  plot = Seurat::DimPlot(data, reduction = reductionType, group.by = annotationField, label = T) +
    patchwork::plot_annotation(titlePlot)
  return(plot)
}

