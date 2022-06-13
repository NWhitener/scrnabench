#' Completes the KMeans Process
#'
#' This function completes the entire kmeans workflow procedure for the log data transformation. This workflow handles data transformation,
#' dimensionality reduction, and clustering. The results will be stored in a csv file called kmeans_clusters.csv
#'
#' @param datasets The names of the datasets for the completion of the workflow
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @param path The file path for saving the results of the workflow
#' @export
complete_kmeans <- function(datasets, seed = 1, path = '.', k = 10)
{
  set.seed(seed)
  file = paste(path, "/kmeans_clusters.csv", sep = '')
  dataList <- extract_datasets(datasets)
  dataList <- preprocess(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_log(dataList)
  dataList <- select_hvg(dataList)
  dataList <- scale_data(dataList)
  dataList <- run_pca(dataList)
  dataList <- run_umap(dataList)
  dataList <- run_tsne(dataList)
  dataList <- run_kmeans(dataList, reductionChoosen = 'pca', k)
  dataList <- run_kmeans(dataList, reductionChoosen = 'umap', k)
  dataList <- run_kmeans(dataList, reductionChoosen = 'tsne', k)
  sil_pca = run_silhouette(dataList, 'pca')
  sil_tsne = run_silhouette(dataList, 'tsne')
  sil_umap = run_silhouette(dataList, 'umap')
  dunn_pca = run_dunn(dataList, 'pca')
  dunn_tsne = run_dunn(dataList, 'tsne')
  dunn_umap = run_dunn(dataList, 'umap')
  ncells = NULL
  nclust = NULL
  for(i in 1:length(names(dataList)))
  {
    ncels_value = length(dataList[[i]]@meta.data$kmeans_cluster_pca)
    nclust_val = length(unique(dataList[[i]]@meta.data$kmeans_cluster_pca))
    ncells = append(ncells,ncels_value)
    nclust = append(nclust, nclust_val)
  }
  results_table = NULL
  results_table = cbind(results_table, sil_pca, sil_umap[,2],sil_tsne[,2],dunn_pca[,2],dunn_umap[,2],dunn_tsne[,2], ncells, nclust )
  colnames(results_table) = c("ID", "Silhouette_PCA", "Silhouette_UMAP", "Silhouette_TSNE", "Dunn_PCA", "Dunn_UMAP", "Dunn_TSNE", "Number of Cells", "Number of Clusters")
  write.table(results_table, file = file, sep = ',', row.names = F, col.names = T)
  return(results_table)
}

#' Completes the Seurat Process
#'
#' This function completes the entire kmeans workflow procedure for the log data transformation. This workflow handles data transformation,
#' dimensionality reduction, and clustering. The results will be stored in a csv file called seurat_clusters.csv
#'
#' @param datasets The names of the datasets for the completion of the workflow
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @param path The file path for saving the results of the workflow
#' @export
complete_seurat <-function(datasets, seed = 1, path = '.', reduction = 'pca')
{
  set.seed(seed)
  file = paste(path, "/seurat_clusters.csv", sep = '')
  dataList <- extract_datasets(datasets)
  print("Extracted")
  dataList <- preprocess(dataList)
  print("Preprocessed")
  dataList <- annotate_datasets(dataList)
  print("Annotated")
  dataList <- run_log(dataList)
  print("Log Normalized")
  dataList <- select_hvg(dataList)
  print("Highly Variable Genes Selected")
  dataList <- scale_data(dataList)
  print("Scaled and Centered")
  dataList <- run_pca(dataList)
  print("Principle Componenet Analysis Completed")
  dataList <- run_umap(dataList)
  print(("Uniform Manifold Approximation and Projection Applied"))
  dataList <- run_tsne(dataList)
  print("T distributed stochastic neighbor embedding")
  dataList <- run_seurat_cluster(dataList, reductionChoosen = reduction)
  print("CLustering done")
  sil_seurat = run_silhouette(dataList, reductionChoosen = reduction, method = "seurat")
  print("SIL calculated")
  dunn_seurat = run_dunn(dataList, reductionChoosen = reduction, method = 'seurat')
  print("DUNN calculated")
  ncells = NULL
  nclust = NULL
  for(i in 1:length(names(dataList)))
  {
    ncels_value = length(dataList[[i]]@meta.data$seurat_cluster)
    nclust_val = length(unique(dataList[[i]]@meta.data$seurat_cluster))
    ncells = append(ncells,ncels_value)
    nclust = append(nclust, nclust_val)
  }
  results_table = NULL
  results_table = cbind(results_table, sil_seurat, dunn_seurat[,2], ncells, nclust)
  colnames(results_table) = c("ID", "Silhouette", "Dunn", "Number of Cells", "Number of Clusters")
  write.table(results_table, file = file, sep = ',' , row.names = F, col.names = T)
  return(results_table)
}

#' Completes the KMeans Process on GFICF data
#'
#' This function completes the entire kmeans workflow procedure for the GFICF data transformation. This workflow handles data transformation,
#' dimensionality reduction, and clustering. The results will be stored in a csv file called kmeans_clusters_gficf.csv
#'
#' @param datasets The names of the datasets for the completion of the workflow
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @param path The file path for saving the results of the workflow
#' @export
complete_kmeans_gficf <- function(datasets, seed = 1, path = '.')
{
  set.seed(seed)
  file = paste(path, "/kmeans_clusters_gficf.csv", sep = '')
  dataList <- extract_datasets(datasets)
  dataList <- run_tfidf(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- select_hvg(dataList)
  dataList <- run_pca(dataList)
  dataList <- run_umap(dataList)
  dataList <- run_tsne(dataList)
  dataList <- run_kmeans(dataList, reductionChoosen = 'pca')
  dataList <- run_kmeans(dataList, reductionChoosen = 'umap')
  dataList <- run_kmeans(dataList, reductionChoosen = 'tsne')
  sil_pca = run_silhouette(dataList, 'pca')
  sil_tsne = run_silhouette(dataList, 'tsne')
  sil_umap = run_silhouette(dataList, 'umap')
  dunn_pca = run_dunn(dataList, 'pca')
  dunn_tsne = run_dunn(dataList, 'tsne')
  dunn_umap = run_dunn(dataList, 'umap')
  ncells = NULL
  nclust = NULL
  for(i in 1:length(names(dataList)))
  {
    ncels_value = length(dataList[[i]]@meta.data$kmeans_cluster_pca)
    nclust_val = length(unique(dataList[[i]]@meta.data$kmeans_cluster_pca))
    ncells = append(ncells,ncels_value)
    nclust = append(nclust, nclust_val)
  }
  results_table = NULL
  results_table = cbind(results_table, sil_pca, sil_umap[,2],sil_tsne[,2],dunn_pca[,2],dunn_umap[,2],dunn_tsne[,2], ncells, nclust )
  colnames(results_table) = c("ID", "Silhouette_PCA", "Silhouette_UMAP", "Silhouette_TSNE", "Dunn_PCA", "Dunn_UMAP", "Dunn_TSNE", "Number of Cells", "Number of Clusters")
  write.table(results_table, file = file, sep = ',', row.names = F, col.names = T)
  return(results_table)
}

#' Completes the Seurat Process on GFICF data
#'
#' This function completes the entire Seurat workflow procedure for the GFICF data transformation. This workflow handles data transformation,
#' dimensionality reduction, and clustering. The results will be stored in a csv file called seurat_clusters_gficf.csv
#'
#' @param datasets The names of the datasets for the completion of the workflow
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @param path The file path for saving the results of the workflow
#' @export
complete_seurat_gficf <- function(datasets, seed = 1, path = '.', reduction = 'pca')
{
  set.seed(seed)
  file = paste(path, "/seurat_clusters_gficf.csv", sep = '')
  dataList <- extract_datasets(datasets)
  dataList <- run_tfidf(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- select_hvg(dataList)
  dataList <- run_pca(dataList)
  dataList <- run_umap(dataList)
  dataList <- run_tsne(dataList)
  dataList <- run_seurat_cluster(dataList, reductionChoosen = reduction)
  sil_seurat = run_silhouette(dataList, reductionChoosen = reduction, method = "seurat")
  dunn_seurat = run_dunn(dataList, reductionChoosen = reduction, method = 'seurat')
  ncells = NULL
  nclust = NULL
  for(i in 1:length(names(dataList)))
  {
    ncels_value = length(dataList[[i]]@meta.data$seurat_cluster)
    nclust_val = length(unique(dataList[[i]]@meta.data$seurat_cluster))
    ncells = append(ncells,ncels_value)
    nclust = append(nclust, nclust_val)
  }
  results_table = NULL
  results_table = cbind(results_table, sil_seurat, dunn_seurat[,2])
  colnames(results_table) = c("ID", "Silhouette", "Dunn", "Number of Cells", "Number of Clusters")
  write.table(results_table, file = file, sep = ',' , row.names = F, col.names = T)
  return(results_table)
}

#' Completes the Harmony Workflow on Log Data
#'
#' This function completes the harmony workflow procedure for the log data transformation for merging datasets. This workflow handles data transformation,
#' dimensionality reduction, merging, and clustering.
#'
#' @param datasets The names of the datasets for the completion of the workflow
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @export
run_harmony_workflow <- function(datasets, seed = 1)
{
  set.seed(seed)
  dataList = extract_datasets(datasets)
  dataList = extract_common_genes(dataList)
  dataList = merge_datasets(dataList)
  dataList = preprocess(dataList)
  dataList = annotate_datasets(dataList)
  dataList = run_log(dataList)
  dataList = select_hvg(dataList)
  dataList = scale_data(dataList)
  dataList = run_pca(dataList)
  dataList = run_harmony(dataList)
  dataList = run_kmeans(dataList, reductionChoosen = 'harmony')
  dataList = run_silhouette(dataList, reductionChoosen = "harmony")
  retunr(dataList)
}
#' Completes the FastMNN Workflow on Log Data
#'
#' This function completes the FastMNN workflow procedure for the log data transformation for merging datasets. This workflow handles data transformation,
#' dimensionality reduction, merging, and clustering.
#'
#' @param datasets The names of the datasets for the completion of the workflow
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @export
run_fastmnn_workflow <- function(datasets, seed = 1)
{
  ##REWRITE
  dataList <- extract_datasets(datasets)
  dataList <- extract_common_genes(dataList)
  dataList <- merge_datasets(dataList)
  dataList <- preprocess(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_log(dataList)
  dataList <- select_hvg(dataList)
  dataList <- run_fastmnn(dataList)
  dataList <- run_umap(dataList)
  dataList <- run_seurat_cluster(dataList)
  return(dataList)
}

#' Completes the CCA Workflow on Log Data
#'
#' This function completes the CCA workflow procedure for the log data transformation for merging datasets. This workflow handles data transformation,
#' dimensionality reduction, merging, and clustering.
#'
#' @param datasets The names of the datasets for the completion of the workflow
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @export
run_cca_workflow <- function(datasets, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(datasets)
  dataList <- preprocess(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_log(dataList)
  dataList <- select_hvg(dataList)
  dataList <- run_cca(dataList)
  dataList <- scale_data(dataList)
  dataList <- run_pca(dataList)
  dataList <- run_umap(dataList)
  dataList <- run_seurat_cluster(dataList)
  return(dataList)
}

#' Completes the SCtransform Workflow on Log Data
#'
#' This function completes the sctransform workflow procedure for the log data transformation for merging datasets. This workflow handles data transformation,
#' dimensionality reduction, merging, and clustering.
#'
#' @param datasets The names of the datasets for the completion of the workflow
#' @param seed The seed of randomization for reproducibility, defaults to 1
#' @export
run_sctransform_workflow <- function(datasets, seed = 1)
{
  set.seed(seed)
  dataList <- extracrt_datasets(dataList)
  dataList <- preprocess(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_sctransform(dataList)
  dataList <- scale_data(dataList)
  dataList <- run_pca(dataList)
  dataList <- run_umap(dataList)
  dataList <- run_seurat_cluster(dataList)
  return(dataList)
}


#' Build Seurat Pipeline with the Columns Permuted
#'
#' This functions builds the Seurat Pipeline, with experimental condition. The initial data set is appended with a the same data list,
#' with permuted columns. This allows the testing of the algorithmic stability of the Seurat Pipeline for the specified data list. The
#' Normalized information Distance is returned. If there is true algorithmic stability this score would be 1
#'
#' @param datasets A data set with genes as rows and cells as columns
#' @return Adjusted Rand Index
build_seurat_columns <- function(datasets, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(datasets)
  dataList2 <- extract_datasets(datasets)
  dataList2 <- permute_columns(dataList2)
  dataList <- append(dataList, dataList2)
  dataList <- preprocess(dataList)
  dataList <- run_log(dataList)
  dataList <- run_seurat(dataList)
  x = merge(dataList[[1]]@meta.data, dataList[[2]]@meta.data, by="row.names", all=TRUE)
  return(aricode::ARI(x$seurat_clusters.x, x$seurat_clusters.y))
}

#' Build Seurat Pipeline with the Rows Permuted
#'
#' This functions builds the Seurat Pipeline, with experimental condition. The initial data set is appended with a the same data list,
#' with permuted rows. This allows the testing of the algorithmic stability of the Seurat Pipeline for the specified data list. The
#' Normalized information Distance is returned. If there is true algorithmic stability this score would be 1
#'
#' @param datasets A data set with genes as rows and cells as columns
#' @return Adjusted Rand Index
build_seurat_rows <- function(datasets, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(datasets)
  dataList2 <- extract_datasets(datasets)
  dataList2 <- permute_rows(dataList2)
  dataList <- append(dataList, dataList2)
  dataList <- preprocess(dataList)
  dataList <- run_log(dataList)
  dataList <- run_seurat(dataList)
  x = merge(dataList[[1]]@meta.data, dataList[[2]]@meta.data, by="row.names", all=TRUE)
  return(aricode::ARI(x$seurat_clusters.x, x$seurat_clusters.y))
}

#' Run the Seurat Pipeline with column Permutations
#'
#' This function runs the built Seurat pipeline with column permutations on a data list. The results are appened into
#' the rows of a new matrix, with the datasets and the ARI recorded
#' @return a results matrix with ARI and datasets number
#' @export
run_seurat_columns <- function()
{
  results = cbind('datasets', 'ARI')
  for (datasets in datasets)
  {
    result <- build_seurat_columns(datasets)
    results <- rbind(results, cbind(datasets, result))
  }
  return(results)
}

#' Run the Seurat Pipeline with row Permutations
#'
#' This function runs the built Seurat pipeline with row permutations on a data list. The results are appended into
#' the rows of a new matrix, with the datasets and the ARI recorded
#' @return a results matrix with ARI and datasets number
#' @export
run_seurat_rows <- function()
{
  results = cbind('datasets', 'ARI')
  for (datasets in datasets)
  {
    result <- build_seurat_rows(datasets)
    results <- rbind(results, cbind(datasets, result))
  }
  return(results)
}


#' Runs the Seurat Pipeline
#'
#' This function completes the Seurat Pipeline on a list of data. This process starts with finding the variable features (genes), then
#' scaling the data, running PCA, UMPA, and then find clusters.  This step can be done to test your benchamrk data sets
#' or for down stream analysis. For more information on the Seurat Pipeline
#' please see https://satijalab.org/seurat/index.html
#'
#' @param dataList A data list of Data sets in the Gene in Row and Cell in Columns format
#' @export
run_seurat <- function(dataList, seed =1)
{
  set.seed(seed)
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- Seurat::RunPCA(x, npcs = 30, verbose = FALSE)
    x <- Seurat::RunUMAP(x, reduction = "pca", dims = 1:30)
    x <- Seurat::FindNeighbors(x, reduction = "pca", dims = 1:30)
    x <- Seurat::FindClusters(x, resolution = 0.5)
  })
  return(dataList)
}



#' Run Integration  (Need a better name)
#'
#' This functions completes an integration pipeline, using
#' 4 integration techniques. It implements the Harmony, FastMNN, CCA, and scTranform integration techniques.
#' These are implemented via the Seurat and Harmony packages. See ?run_harmony, ?run_cca, ?run_fastmnn, and ?run_sctransform
#' for more information on these pipeline structures. Each of these integration will write the integrated data set to a file
#' in the current directory
#'
#' @param datasets A list of data set names that you would like to integrate with
#' the 4 pipelines
#' @export
run_integration <- function(datasets)
{
  # Harmony
  data <- run_harmony(datasets, batch_name = "ID")
  #write_output(data[1], 'harmony')

 #FastMnn
  data <- run_fastmnn(datasets, batch_column)
  #write_output(data, 'fastmnn')

  # CCA
  data <- run_cca(datasets)
  #write_output(data, 'cca')

  # ScTransform
  data <- run_sctransform(datasets)
  #write_output(data, 'sctransform')
}


#' Run Integration Duplicated (Need a better name)
#'
#' This functions completes an integration pipeline on a dupilcated set of the original data, using
#' 4 integration techniques. It implements the Harmony, FastMNN, CCA, and scTranform integration techniques.
#' These are implemented via the Seurat and Harmony packages. See ?run_harmony, ?run_cca, ?run_fastmnn, and ?run_sctransform
#' for more information on these pipeline structures. Each of these integration will write the integrated data set to a file
#' in the current directory
#'
#' @param datasets A list of data set names that you would like to integrate with
#' the 4 pipelines
#' @param dup The number of duplications that would like to execute
#' @export
run_duplicate_integrations <- function(datasets, dup)
{
  # Harmony
  dataList <- duplicate_datasets(datasets, dup)
  data <- run_harmony(data, batch_column)
  write_output(data, 'harmony')

  # FastMNN
  dataList <- duplicate_datasets(datasets, dup)
  data <- run_fastmnn(data, batch_column)
  write_output(data, 'fastmnn')

  # CCA
  dataList <- duplicate_datasets(datasets, dup)
  data <- run_cca(dataList)
  write_output(data, 'cca')

  # ScTransform
  dataList <- duplicate_datasets(datasets, dup)
  data <- run_sctransform(dataList)

  write_output(data, 'sctransform')
}

#' Write Output
#'
#' This function writes the output of the data integration pipelines
#' to a file with the associated labels, and data information.
#'
#' @param data The data the needs to be output
#' @param prefix The prefix of the integratio technique used
write_output <- function(data, prefix)
{
  print(prefix)
  print(dim(as.data.frame(Seurat::GetAssayData(data))[Seurat::VariableFeatures(data),]))
  path = paste("./", prefix, sep="")

  write.csv(t(as.data.frame(Seurat::GetAssayData(data))[Seurat::VariableFeatures(data),]), paste(path, "_hvg.csv", sep=""), quote=F)
  write.csv(Seurat::Embeddings(data, 'umap'), paste(path, "_umap.csv", sep=""), quote=F)
  if(prefix=='harmony')
  {
    write.csv(Seurat::Embeddings(data, 'harmony'), paste(path, "_harmonized.csv", sep=""), quote=F)
  }
  else if (prefix == 'fastmnn')
  {
    write.csv(Seurat::Embeddings(data, 'mnn'), paste(path, "_harmonized.csv", sep=""), quote=F)
  }
  else
  {
    write.csv(Seurat::Embeddings(data, 'pca'), paste(path, "_pca.csv", sep=""), quote=F)
  }
  write.csv(data@meta.data, paste(path, "_labels.csv", sep=""), quote=F)
}




#' Get Scenario
#'
#' This function gets the Scenario used in the creation of the data and the sequencing
#' technique used in collection. Only to be used internally.
#'
#' @param scenario_id The scenario_id that you would like to identify
get_scenario <- function(scenario_id)
{
  scenarios = list(

    # CELL_LINE classification scenario 1 to 21
    c('10X_LLU_A_cellranger2.0', '10X_LLU_B_cellranger2.0'),
    c('10X_LLU_A_cellranger3.1', '10X_LLU_B_cellranger3.1'),
    c('10X_LLU_A_zumi', '10X_LLU_B_zumi'),
    c('10X_LLU_A_umitools', '10X_LLU_B_umitools'),

    c('10X_NCI_A_cellranger2.0', '10X_NCI_B_cellranger2.0'),
    c('10X_NCI_A_cellranger3.1', '10X_NCI_B_cellranger3.1'),
    c('10X_NCI_A_zumi', '10X_NCI_B_zumi'),
    c('10X_NCI_A_umitools', '10X_NCI_B_umitools'),

    c('10X_NCI_M_A_cellranger2.0', '10X_NCI_M_B_cellranger2.0'),
    c('10X_NCI_M_A_cellranger3.1', '10X_NCI_M_B_cellranger3.1'),
    c('10X_NCI_M_A_zumi', '10X_NCI_M_B_zumi'),
    c('10X_NCI_M_A_umitools', '10X_NCI_M_B_umitools'),

    c('C1_FDA_HT_A_featureCounts', 'C1_FDA_HT_B_featureCounts'),
    c('C1_FDA_HT_A_rsem', 'C1_FDA_HT_B_rsem'),
    c('C1_FDA_HT_A_kallisto', 'C1_FDA_HT_B_kallisto'),

    c('ICELL8_PE_A_featureCounts', 'ICELL8_PE_B_featureCounts'),
    c('ICELL8_PE_A_rsem', 'ICELL8_PE_B_rsem'),
    c('ICELL8_PE_A_kallisto', 'ICELL8_PE_B_kallisto'),

    c('ICELL8_SE_A_featureCounts', 'ICELL8_SE_B_featureCounts'),
    c('ICELL8_SE_A_rsem', 'ICELL8_SE_B_rsem'),
    c('ICELL8_SE_A_kallisto', 'ICELL8_SE_B_kallisto'),

    # CELL_LINE classification 22 - 24
    c('C1_LLU_A_featureCounts', 'C1_LLU_B_featureCounts'),
    c('C1_LLU_A_rsem', 'C1_LLU_B_rsem'),
    c('C1_LLU_A_kallisto', 'C1_LLU_B_kallisto'),

    # Bias detection scenario 25 to 27
    c('10X_LLU_A_cellranger2.0', '10X_LLU_A_cellranger3.1', '10X_LLU_A_zumi', '10X_LLU_A_umitools',
      '10X_NCI_A_cellranger2.0', '10X_NCI_A_cellranger3.1', '10X_NCI_A_zumi', '10X_NCI_A_umitools',
      '10X_NCI_M_A_cellranger2.0', '10X_NCI_M_A_cellranger3.1', '10X_NCI_M_A_zumi',
      '10X_NCI_M_A_umitools',
      'ICELL8_PE_A_featureCounts',
      'ICELL8_PE_A_rsem',
      'ICELL8_PE_A_kallisto',
      'ICELL8_SE_A_featureCounts',
      'ICELL8_SE_A_rsem',
      'ICELL8_SE_A_kallisto',
      'C1_FDA_HT_A_featureCounts',
      'C1_FDA_HT_A_rsem',
      'C1_FDA_HT_A_kallisto'),

    c('10X_LLU_B_cellranger2.0', '10X_LLU_B_cellranger3.1', '10X_LLU_B_zumi', '10X_LLU_B_umitools',
      '10X_NCI_B_cellranger2.0', '10X_NCI_B_cellranger3.1', '10X_NCI_B_zumi', '10X_NCI_B_umitools',
      '10X_NCI_M_B_cellranger2.0', '10X_NCI_M_B_cellranger3.1',
      '10X_NCI_M_B_zumi', '10X_NCI_M_B_umitools',
      'ICELL8_PE_B_featureCounts',
      'ICELL8_PE_B_rsem',
      'ICELL8_PE_B_kallisto',
      'ICELL8_SE_B_featureCounts',
      'ICELL8_SE_B_rsem',
      'ICELL8_SE_B_kallisto',
      'C1_FDA_HT_B_featureCounts',
      'C1_FDA_HT_B_rsem',
      'C1_FDA_HT_B_kallisto'),

    c('10X_LLU_A_cellranger2.0', '10X_LLU_A_cellranger3.1', '10X_LLU_A_zumi', '10X_LLU_A_umitools',
      '10X_LLU_B_cellranger2.0', '10X_LLU_B_cellranger3.1', '10X_LLU_B_zumi', '10X_LLU_B_umitools',
      '10X_NCI_A_cellranger2.0', '10X_NCI_A_cellranger3.1', '10X_NCI_A_zumi', '10X_NCI_A_umitools',
      '10X_NCI_B_cellranger2.0', '10X_NCI_B_cellranger3.1', '10X_NCI_B_zumi', '10X_NCI_B_umitools',
      '10X_NCI_M_A_cellranger2.0', '10X_NCI_M_A_cellranger3.1', '10X_NCI_M_A_zumi',
      '10X_NCI_M_A_umitools', '10X_NCI_M_B_cellranger2.0', '10X_NCI_M_B_cellranger3.1',
      '10X_NCI_M_B_zumi', '10X_NCI_M_B_umitools',
      'ICELL8_PE_A_featureCounts', 'ICELL8_PE_B_featureCounts',
      'ICELL8_PE_A_rsem', 'ICELL8_PE_B_rsem',
      'ICELL8_PE_A_kallisto', 'ICELL8_PE_B_kallisto',
      'ICELL8_SE_A_featureCounts', 'ICELL8_SE_B_featureCounts',
      'ICELL8_SE_A_rsem', 'ICELL8_SE_B_rsem',
      'ICELL8_SE_A_kallisto', 'ICELL8_SE_B_kallisto',
      'C1_FDA_HT_A_featureCounts', 'C1_FDA_HT_B_featureCounts',
      'C1_FDA_HT_A_rsem', 'C1_FDA_HT_B_rsem',
      'C1_FDA_HT_A_kallisto', 'C1_FDA_HT_B_kallisto'),

    # Bias detection scenario 28 to 48
    c('10X_LLU_A_cellranger2.0', '10X_LLU_B_cellranger2.0'),
    c('10X_LLU_A_cellranger3.1', '10X_LLU_B_cellranger3.1'),
    c('10X_LLU_A_zumi', '10X_LLU_B_zumi'),
    c('10X_LLU_A_umitools', '10X_LLU_B_umitools'),

    c('10X_NCI_A_cellranger2.0', '10X_NCI_B_cellranger2.0'),
    c('10X_NCI_A_cellranger3.1', '10X_NCI_B_cellranger3.1'),
    c('10X_NCI_A_zumi', '10X_NCI_B_zumi'),
    c('10X_NCI_A_umitools', '10X_NCI_B_umitools'),

    c('10X_NCI_M_A_cellranger2.0', '10X_NCI_M_B_cellranger2.0'),
    c('10X_NCI_M_A_cellranger3.1', '10X_NCI_M_B_cellranger3.1'),
    c('10X_NCI_M_A_zumi', '10X_NCI_M_B_zumi'),
    c('10X_NCI_M_A_umitools', '10X_NCI_M_B_umitools'),

    c('ICELL8_PE_A_featureCounts', 'ICELL8_PE_B_featureCounts'),
    c('ICELL8_PE_A_rsem', 'ICELL8_PE_B_rsem'),
    c('ICELL8_PE_A_kallisto', 'ICELL8_PE_B_kallisto'),

    c('ICELL8_SE_A_featureCounts', 'ICELL8_SE_B_featureCounts'),
    c('ICELL8_SE_A_rsem', 'ICELL8_SE_B_rsem'),
    c('ICELL8_SE_A_kallisto', 'ICELL8_SE_B_kallisto'),

    c('C1_FDA_HT_A_featureCounts', 'C1_FDA_HT_B_featureCounts'),
    c('C1_FDA_HT_A_rsem', 'C1_FDA_HT_B_rsem'),
    c('C1_FDA_HT_A_kallisto', 'C1_FDA_HT_B_kallisto'),

    # 49 - 69
    c('10X_LLU_A_cellranger2.0', '10X_LLU_A_cellranger2.0'),
    c('10X_LLU_A_cellranger3.1', '10X_LLU_A_cellranger3.1'),
    c('10X_LLU_A_zumi', '10X_LLU_A_zumi'),
    c('10X_LLU_A_umitools', '10X_LLU_A_umitools'),

    c('10X_NCI_A_cellranger2.0', '10X_NCI_A_cellranger2.0'),
    c('10X_NCI_A_cellranger3.1', '10X_NCI_A_cellranger3.1'),
    c('10X_NCI_A_zumi', '10X_NCI_A_zumi'),
    c('10X_NCI_A_umitools', '10X_NCI_A_umitools'),

    c('10X_NCI_M_A_cellranger2.0', '10X_NCI_M_A_cellranger2.0'),
    c('10X_NCI_M_A_cellranger3.1', '10X_NCI_M_A_cellranger3.1'),
    c('10X_NCI_M_A_zumi', '10X_NCI_M_A_zumi'),
    c('10X_NCI_M_A_umitools', '10X_NCI_M_A_umitools'),

    c('ICELL8_PE_A_featureCounts', 'ICELL8_PE_A_featureCounts'),
    c('ICELL8_PE_A_rsem', 'ICELL8_PE_A_rsem'),
    c('ICELL8_PE_A_kallisto', 'ICELL8_PE_A_kallisto'),

    c('ICELL8_SE_A_featureCounts', 'ICELL8_SE_A_featureCounts'),
    c('ICELL8_SE_A_rsem', 'ICELL8_SE_A_rsem'),
    c('ICELL8_SE_A_kallisto', 'ICELL8_SE_A_kallisto'),

    c('C1_FDA_HT_A_featureCounts', 'C1_FDA_HT_A_featureCounts'),
    c('C1_FDA_HT_A_rsem', 'C1_FDA_HT_A_rsem'),
    c('C1_FDA_HT_A_kallisto', 'C1_FDA_HT_A_kallisto')
  )

  return(scenarios[[scenario_id]])
}




