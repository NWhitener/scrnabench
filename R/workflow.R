#' Completes the clustering workflow
#'
#' This function completes the entire clustering workflow, on different data transformation types,
#' and with different clustering methods. The log data transformation pipeline includes filter_data, annotating,
#' applying the log transformation, selecting the Highly Variable Genes, and finally scaling the data before
#' dimensionality reduction, and clustering. The TFIDF data transformation completes the TFIDF data transformation,
#' annotates the data, and selects the Highly Variable Genes before applying dimensionality reduction and finally
#' clustering. The clustering parameters are set to default, except for the number of clusters in
#' kmeans clustering, which is set to 10
#'
#' @param dataList A list of the data sets
#' @param method The method of clustering to use, defaults to Kmeans
#' @param transformationType The data transformation to be used, defaults to the log transformation
#' @param seed The seed to be set for reproducibility, defaults to 1
#' @param numberClusters The number of clusters to use for Kmeans clustering, defaults to 10
#' @return the Data list with all components of the workflow completed
#' @export
run_clustering_workflow <- function(dataList, method = 'kmeans', transformationType = 'log', seed = 1, numberClusters = 10)
{
  if(is.list(dataList))
  {
  set.seed(seed)

  if(transformationType == 'log')
  {
    dataList <- filter_data(dataList)
    dataList <- annotate_datasets(dataList)
    dataList <- run_log(dataList)
    dataList <- select_hvg(dataList)
    dataList <- scale_data(dataList)
  }
  else if (transformationType == 'tfidf')
  {

    dataList <- filter_data(dataList)
    dataList <- run_tfidf(dataList)
    dataList <- annotate_datasets(dataList)
    dataList <- select_hvg(dataList)
  }
  else
  {
    stop("Unknown transformation type. Only log and tfidf transformations are supported.")
  }

  dataList <- run_pca(dataList, numComponents = 10)
  dataList <- run_umap(dataList, numDimensions = 10)
  dataList <- run_tsne(dataList, numDimensions = 10)

  if(method == 'kmeans')
  {
    dataList <- run_kmeans_clustering(dataList, reductionType = 'pca', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'umap', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'tsne', numberClusters)
  }
  else if(method == 'seurat')
  {
    dataList <- run_seurat_clustering(dataList, reductionType = 'pca', numberComponents = 10)
    dataList <- run_seurat_clustering(dataList, reductionType = 'umap', numberComponents = 2)
    dataList <- run_seurat_clustering(dataList, reductionType = 'tsne', numberComponents = 2)
  }
  else
  {
    stop("Unknown clustering method. Only seurat and kmeans methods are supported.")
  }

  print(create_internal_cluster_validation_report(dataList, method = method))
  }
  else
  {
    stop("A data list of datasets is required to run the clustering workflow on the datasets")
  }
  return(dataList)

}

#' Completes the Harmony Workflow
#'
#' This function completes the harmony workflow procedure for the log data transformation for merging datasets. This workflow handles data transformation,
#' dimensionality reduction, merging, and clustering.
#'
#' @param dataList A list of data sets to merge
#' @param method The clustering method to use
#' @param numberClusters The number of clusters to use in the kmeans clustering
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @return A list of datasets with the harmony workflow applied
#' @export
run_harmony_integration_workflow <- function(dataList, method = 'kmeans', numberClusters = 10, seed = 1)
{
  if(is.list(dataList))
  {
  set.seed(seed)
  dataList = extract_common_genes(dataList)
  dataList = merge_datasets(dataList)
  dataList = filter_data(dataList)
  dataList = annotate_datasets(dataList)
  dataList = run_log(dataList)
  dataList = select_hvg(dataList)
  dataList = scale_data(dataList)
  dataList = run_pca(dataList, numComponents = 10)
  dataList = run_harmony(dataList)
  dataList <- run_umap(dataList, reductionType = 'harmony', numDimensions = 10)
  dataList <- run_tsne(dataList, reductionType = 'harmony', numDimensions = 10)

  if(method == 'kmeans')
  {
    dataList <- run_kmeans_clustering(dataList, reductionType = 'pca', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'harmony', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'umap', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'tsne', numberClusters)
  }
  else if(method == 'seurat')
  {
    dataList <- run_seurat_clustering(dataList, reductionType = 'pca', numberComponents = 10)
    dataList <- run_seurat_clustering(dataList, reductionType = 'harmony', numberComponents = 10)
    dataList <- run_seurat_clustering(dataList, reductionType = 'umap', numberComponents = 2)
    dataList <- run_seurat_clustering(dataList, reductionType = 'tsne', numberComponents = 2)
  }
  else
  {
    stop("Unknown clustering method. Only seurat and kmeans methods are supported.")
  }

  print(create_internal_cluster_validation_report(dataList, method = method))
  }
  else
  {
    stop("A data list of datasets is required to apply the harmony workflow to the datasets")
  }
  return(dataList)
}

#' Completes the FastMNN Workflow
#'
#' This function completes the FastMNN workflow procedure for the log data transformation for merging datasets. This workflow handles data transformation,
#' dimensionality reduction, merging, and clustering.
#'
#' @param dataList A list of data sets to merge
#' @param method The clustering method to use
#' @param numberClusters The number of clusters to use in the kmeans clustering
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @return A list of datasets with the harmony workflow applied
#' @export
run_fastmnn_integration_workflow <- function(dataList, method = 'kmeans', numberClusters = 10, seed = 1)
{
  if(is.list(dataList))
  {
  dataList <- extract_common_genes(dataList)
  dataList <- merge_datasets(dataList)
  dataList <- filter_data(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_log(dataList)
  dataList <- select_hvg(dataList)
  dataList <- run_fastmnn(dataList)
  dataList <- run_umap(dataList, reductionType = 'mnn', numDimensions = 10)
  dataList <- run_tsne(dataList, reductionType = 'mnn', numDimensions = 10)

  if(method == 'kmeans')
  {
    dataList <- run_kmeans_clustering(dataList, reductionType = 'mnn', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'umap', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'tsne', numberClusters)
  }
  else if(method == 'seurat')
  {
    dataList <- run_seurat_clustering(dataList, reductionType = 'mnn', numberComponents = 10)
    dataList <- run_seurat_clustering(dataList, reductionType = 'umap', numberComponents = 2)
    dataList <- run_seurat_clustering(dataList, reductionType = 'tsne', numberComponents = 2)
  }
  else
  {
    stop("Unknown clustering method. Only seurat and kmeans methods are supported.")
  }

  print(create_internal_cluster_validation_report(dataList, method = method))
}
  else
  {
    stop("A data list of datasets is required to apply the fastmnn workflow to the datasets")
  }
  return(dataList)
}

#' Completes the CCA Workflow
#'
#' This function completes the CCA workflow procedure for the log data transformation for merging datasets. This workflow handles data transformation,
#' dimensionality reduction, merging, and clustering.
#'
#' @param dataList A list of data sets to merge
#' @param method The clustering method to use
#' @param numberClusters The number of clusters to use in the kmeans clustering
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @return A list of datasets with the harmony workflow applied
#' @export
run_cca_integration_workflow <- function(dataList, method = 'kmeans', numberClusters = 10, seed = 1)
{
  if(is.list(dataList))
  {
  set.seed(seed)
  dataList <- filter_data(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_log(dataList)
  dataList <- select_hvg(dataList)
  dataList <- run_cca(dataList)
  dataList <- scale_data(dataList)
  dataList <- run_pca(dataList, numComponents = 10)
  dataList <- run_umap(dataList, reductionType = 'pca', numDimensions = 10)
  dataList <- run_tsne(dataList, reductionType = 'pca', numDimensions = 10)

  if(method == 'kmeans')
  {
    dataList <- run_kmeans_clustering(dataList, reductionType = 'pca', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'umap', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'tsne', numberClusters)
  }
  else if(method == 'seurat')
  {
    dataList <- run_seurat_clustering(dataList, reductionType = 'pca', numberComponents = 10)
    dataList <- run_seurat_clustering(dataList, reductionType = 'umap', numberComponents = 2)
    dataList <- run_seurat_clustering(dataList, reductionType = 'tsne', numberComponents = 2)
  }
  else
  {
    stop("Unknown clustering method. Only seurat and kmeans methods are supported.")
  }

  print(create_internal_cluster_validation_report(dataList, method = method))
  }
  else
  {
    stop("A data list of datasets is required to apply the cca workflow to the datasets")
  }
  return(dataList)
}

#' Completes the SCtransform Workflow
#'
#' This function completes the SCtansform workflow procedure for the log data transformation for merging datasets. This workflow handles data transformation,
#' dimensionality reduction, merging, and clustering.
#'
#' @param dataList A list of data sets to merge
#' @param method The clustering method to use
#' @param numberClusters The number of clusters to use in the kmeans clustering
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @return A list of datasets with the harmony workflow applied
#' @export
run_sctransform_integration_workflow <- function(dataList, method = 'kmeans', numberClusters = 10, seed = 1)
{
  if(is.list(dataList))
  {
  set.seed(seed)
  dataList <- filter_data(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_sctransform(dataList)
  dataList <- scale_data(dataList)

  dataList <- run_pca(dataList, numComponents = 10)
  dataList <- run_umap(dataList, reductionType = 'pca', numDimensions = 10)
  dataList <- run_tsne(dataList, reductionType = 'pca', numDimensions = 10)

  if(method == 'kmeans')
  {
    dataList <- run_kmeans_clustering(dataList, reductionType = 'pca', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'umap', numberClusters)
    dataList <- run_kmeans_clustering(dataList, reductionType = 'tsne', numberClusters)
  }
  else if(method == 'seurat')
  {
    dataList <- run_seurat_clustering(dataList, reductionType = 'pca', numberComponents = 10)
    dataList <- run_seurat_clustering(dataList, reductionType = 'umap', numberComponents = 2)
    dataList <- run_seurat_clustering(dataList, reductionType = 'tsne', numberComponents = 2)
  }
  else
  {
    stop("Unknown clustering method. Only seurat and kmeans methods are supported.")
  }

  print(create_internal_cluster_validation_report(dataList, method = method))
  }
  else
  {
    stop("A data list of datasets is required to apply the cca workflow to the datasets")
  }
  return(dataList)
}


#' Completes the Metamorphic Testing workflow
#'
#' This function completes the entire metamorphic testing clustering workflow, on different data transformation types,
#' and with different clustering methods. This tests different types of metamorphic tests, including permutations, adding genes, etc. The original data
#' is clustered using the clustering workflow, and then the each of the metamorphic tests is applied and the clustering is recomputed. The log data transformation pipeline includes filter_data, annotating,
#' applying the log transformation, selecting the Highly Variable Genes, and finally scaling the data before
#' dimensionality reduction, and clustering. The TFIDF data transformation completes the TFIDF data transformation,
#' annotates the data, and selects the Highly Variable Genes before applying dimensionality reduction and finally
#' clustering. The clustering parameters are set to default, except for the number of clusters in
#' kmeans clustering, which is set to 10. If a score of 0 is returned it means that the dataset did not pass any of the metamorphic
#' tests, if  a 1 is returned then it means that the datasets passed the metamorphic test. Any number inbetween 0 and 1 represnts the percent of
#' reductions that passed the metamorphic test
#'
#' @param dataList A list of the data sets
#' @param method The method of clustering to use, defaults to Kmeans
#' @param transformationType The data transformation to be used, defaults to the log transformation
#' @param seed The seed to be set for reproducibility, defaults to 1
#' @param numberClusters The number of clusters to use for Kmeans clustering, defaults to 10
#' @param metamorphicTests a vector of the metamorphic tests to apply, the options are Permute Cells, Modify Gene Counts, Add Duplicate Cell, Permute Genes, Add Zero Variance Gene, Flip Gene Counts
#' @return A report of which metamorphic tests were passed.
#' @export
run_metamorphic_test_worflow <-function(dataList, method = 'kmeans', transformationType = 'log',  seed = 1, numberClusters = 10, metamorphicTests = c(1:6))
{
  metamorphicTestingReport <- NULL
  if(is.list(dataList))
  {
    set.seed(seed)
    perturbations <- c('Permute Cells', 'Modify Gene Counts', 'Add Duplicate Cell', 'Permute Genes', 'Add Zero Variance Gene', 'Flip Gene Counts')
    metamorphicReportList <- NULL
    i <- 1
    metamorphicReportList <- NULL
    metamorphicReportList[[i]] <- create_internal_cluster_validation_report(run_clustering_workflow(dataList,
                                                                                                    method, transformationType,seed, numberClusters), method)
    for(test in metamorphicTests)
    {
      i <- i + 1
      metamorphicReportList[[i]] <- run_metamorphic_test(dataList, test, method, transformationType, seed, numberClusters)
    }
    names(metamorphicReportList) <- c('original', perturbations[metamorphicTests])

    metamorphicTestingReport <- create_metamorphic_testing_report(metamorphicReportList)
  }
  else
  {
    stop("A data list of datasets is required to run metamorphic tests")
  }
  cat(" \n")
  cat("\n Metamorphic Testing Report: \n")
  cat(" \n")
  print(metamorphicTestingReport)
  return(metamorphicTestingReport)
}

