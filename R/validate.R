run_silhoutte <- function(dataList, reduction_choosen = 'pca')
{
  list_return = NULL
  for (i in range(1:length(names(dataList))))
  {
      c = paste("kmeans_cluster_", reduction_choosen, sep ="")
      x = dataList[[i]][[c]]
      y = Seurat::Embeddings(dataList[[i]], reduction = reduction_choosen)
      z <- cluster::silhouette(x[,1],dist(y, "euclidean"))
      x = summary(z)
      list_return = append(list_return, x$avg.width)
      print(x$avg.width)
  }

}


run_dunn <- function(dataList){

  list_return = NULL
  for (i in range(1:length(names(dataList))))
  {
    c = paste("kmeans_cluster_", reduction_choosen, sep ="")
    y = dataList[[i]][[c]]
    x = Seurat::Embeddings(dataList[[i]], reduction = 'pca')
    dunn = clValid::dunn(clusters = y[,1], Data = x)
    list_return = append(list_return, dunn)
    print(dunn)
    }
    return(list_return)
  }

create_kmeans_plot <- function(data, reduction_choosen = 'pca'){
    c = paste("kmeans_cluster_", reduction_choosen, sep ="")
    return(Seurat::DimPlot(data, reduction =reduction_choosen, group.by = c))
}
