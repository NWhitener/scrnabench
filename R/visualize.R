#' Visualization method for the Kmeans Clustering
#'
#' Produces a plot of the Kmeans clustering, based on the dimensionality reduction method choosen
#'
#' @param data An individual data object that you want to visualize
#' @param reduction_choosen The dimensionality reduction method choosen
#' @export
create_kmeans_plot <- function(data, reduction_choosen = 'pca'){
  c = paste("kmeans_cluster_", reduction_choosen, sep ="")
  return(Seurat::DimPlot(data, reduction =reduction_choosen, group.by = c))
}


#' Visualization method for the Seurat Clustering
#'
#' Produces a plot of the Seurat clustering, based on the dimensionality reduction method choosen
#'
#' @param data An individual data object that you want to visualize
#' @param reduction_choosen The dimensionality reduction method choosen
#' @export
create_seurat_plot <- function(data, reduction_choosen = 'pca'){
  return(Seurat::DimPlot(data, reduction =reduction_choosen, group.by = 'seurat_clusters'))
}
