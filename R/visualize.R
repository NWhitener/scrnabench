#' Visualization method for the Clustering Methods
#'
#' Produces a plot of the clustering, based on the dimensionality reduction method chosen. The plot is a Seurat DimPlot
#'
#' @param data An individual data object that you want to visualize
#' @param reductionChoosen The dimensionality reduction method choosen
#' @param method The cluster method to use
#' @export
plot_clusters <- function(data, reductionChoosen = 'pca', method = 'kmeans'){
  if(method == 'kmeans')
  {
  reductionType = paste('kmeans_cluster_', reductionChoosen, sep = '')
  }
  else {
    reductionType = 'seurat_clusters'
  }
  return(Seurat::DimPlot(data, reduction = reductionChoosen, group.by = reductionType))
}

