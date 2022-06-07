#' COmpletes the KMeans Process
#'
#' @export
complete_kmeans <- function(datasets, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(datasets)
  #dataList <- extract_common_genes(dataList)
  #dataList <- merge_datasets(dataList)
  dataList <- preprocess(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_log(dataList)
  dataList <- select_hvg(dataList)
  dataList <- scale_data(dataList)
  dataList <- run_pca(dataList)
  dataList <- run_umap(dataList)
  dataList <- run_tsne(dataList)
  dataList <- run_kmeans(dataList, reduction_choosen = 'pca')
  dataList <- run_kmeans(dataList, reduction_choosen = 'umap')
  dataList <- run_kmeans(dataList, reduction_choosen = 'tsne')
  x = run_silhouette(dataList, 'pca')
  y = run_silhouette(dataList, 'tsne')
  z = run_silhouette(dataList, 'umap')
  w = run_dunn(dataList, 'pca')
  r = run_dunn(dataList, 'tsne')
  f = run_dunn(dataList, 'umap')
  results_table = NULL
  results_table = cbind(results_table, x, y[,2],z[,2],w[,2],r[,2],f[,2])
  colnames(results_table) = c("ID", "Silhouette_PCA", "Silhouette_UMAP", "Silhouette_TSNE", "Dunn_PCA", "Dunn_UMAP", "Dunn_TSNE")
  write.table(results_table, file = "/Users/nathanwhitener/kmeans.csv", sep = ',')
  return(results_table)
}


complete_seurat <-function(datasets, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(datasets)
  #dataList <- extract_common_genes(dataList)
  #dataList <- merge_datasets(dataList)
  dataList <- preprocess(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_log(dataList)
  dataList <- select_hvg(dataList)
  dataList <- scale_data(dataList)
  dataList <- run_pca(dataList)
  dataList <- run_umap(dataList)
  dataList <- run_tsne(dataList)
  dataList <- run_cluster(dataList)
  x = run_silhouette(dataList, reduction_choosen = "pca", method = "seurat")
  y = run_dunn(dataList, reduction_choosen = 'pca', method = 'seurat')
  results_table = NULL
  results_table = cbind(results_table, x, y[,2])
  print("Here")
  print(results_table)
  colnames(results_table) = c("ID", "Silhouette_PCA", "Dunn_PCA")
  write.table(results_table, file = "/Users/nathanwhitener/seurat.csv", sep = ',')
  return(results_table)
}

#' Run CCA
#'
#' This function runs the data integration protocol detailed in the Seurat "Introduction to scRNA-seq
#' integration" found at https://satijalab.org/seurat/articles/integration_introduction.html. This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param dataList A data list of data sets to integrate using the cca protocol
#' @return data.combined A data list of the combined data from the cca protocol
#' @export
run_cca <- function(idx,seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(idx)
  dataList <- preprocess(dataList)
  #dataList <- run_gficf(dataList)
  dataList <- run_log(dataList)

  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- select_hvg(x)
  })

  features <- Seurat::SelectIntegrationFeatures(object.list = dataList)
  k.filter <- min(200, min(sapply(dataList, ncol)))
  data.anchors <- Seurat::FindIntegrationAnchors(object.list = dataList, anchor.features = features, k.filter=k.filter)
  data.combined <- Seurat::IntegrateData(anchorset = data.anchors)
  SeuratObject::DefaultAssay(data.combined) <- "integrated"

  data.combined <- scale_data(data.combined)
  data.combined <- run_pca(data.combined)
  data.combined <- run_umap(data.combnined)
  data.combined <- run_cluster(data.combined)

  return(data.combined)
}

#' Run Harmony
#'
#' This function runs the data integration protocol detailed in the Seurat ## FIND VING. AND LINK . This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param dataList A data list of data sets to integrate using the harmony protocol
#' @return data.combined A data list of the combined data from the harmony protocol
#' @export
run_harmony <- function(dataList, batch_name='SID', seed = 1)
{
  ##STOP ON LESS THAN 2 DATASETS
  set.seed(seed)
  if(length(dataList) < 2 )
  {
    stop("Only two dataset provided, more are required")
    return(dataList)
  }
  dataList <- lapply(X = dataList, FUN = function(x) {
    x <- harmony::RunHarmony(object = x, group.by.vars = batch_name)
  })
  return(dataList)
}

#' Run fastmnn
#'
#' This function runs the data integration protocol detailed in the Seurat ## FIND VING. AND LINK . This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param dataList A data list of data sets to integrate using the fastmnn protocol
#' @return data.combined A data list of the combined data from the fastmnn protocol
#' @export
run_fastmnn <- function(idx, batch_name, seed =1)
{

  ##REWRITE
  set.seed(seed)
  dataList <- extract_datasets(idx)
  dataList <- extract_common_genes(dataList)
  dataList <- merge_datasets(dataList)
  dataList <- preprocess(dataList)
  dataList <- run_log(dataList) #LOG
  #dataList<- run_gficf(dataList)
  dataList <- select_hvg(dataList)
  dataList <- SeuratWrapper::RunFastMNN(object.list = Seurat::SplitObject(dataList, split.by = batch_name))
  dataList <- run_umap(dataList)
  dataList <- run_cluster(dataList)
  SeuratObject::DefaultAssay(dataList) <- "mnn.reconstructed"
  return(dataList)

}


#' Run sctransform
#'
#' This function runs the data integration protocol detailed in the Seurat "Using sctransform in Seurtat"
#' found at https://satijalab.org/seurat/articles/sctransform_vignette.html . This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param dataList A data list of data sets to integrate using the sctransform protocol
#' @return data.combined A data list of the combined data from the sctransform protocol
#' @export
run_sctransform <- function(dataList, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(idx)
  dataList <- preprocess(dataList)
  dataList <- lapply(X = dataList, FUN = Seurat::SCTransform)
  features <- Seurat::SelectIntegrationFeatures(object.list = dataList, nfeatures=2000)
  dataList <- Seurat::PrepSCTIntegration(object.list = dataList, anchor.features = features)
  k.filter <- min(200, min(sapply(dataList, ncol)))
  data.anchors <- Seurat::FindIntegrationAnchors(object.list = dataList, normalization.method = "SCT", anchor.features = features, k.filter = k.filter)
  data.combined <- Seurat::IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
  SeuratObject::DefaultAssay(data.combined) <- "integrated"
  data.combined <- scale_data(data.combined)
  data.combined <- run_pca(data.combined)
  data.combined <- run_umap(data.combined)
  data.combined <- run_cluster(data.combined)
  return(data.combined)

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

#' Log Workflow Seurat
#'
#' This workflow is designed to apply the log data transformation on the datasets
#' and test it's ability to handle pca, umap, and seurat clustering.
#'
#' @param idx The names of the datasets
#' @param seed The random seed that you would like to use on the dataset
#' @return a result list with the dataset ID and the cluster it belongs to
#' @export
log_workflow_seurat <- function(idx, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(idx)
  #dataList <- extract_common_genes(dataList)
  dataList <- preprocess(dataList)
  dataList <- annotate_datasets(dataList)
  dataList <- run_log(dataList)
  dataList <- select_hvg(dataList)
  dataList <- scale_data(dataList)
  dataList <- run_pca(dataList)
  dataList <- run_umap(dataList)
  dataList <- run_cluster(dataList)
  result <- cbind(dataList[[1]]@meta.data$ID[1],
                  length(dataList[[1]]@meta.data$seurat_clusters),
                  length(unique(dataList[[1]]@meta.data$seurat_clusters)))
  return(result)
}


#' GFICF Workflow Seurat
#'
#' This workflow is designed to apply the GFICF data transformation on the datasets
#' and test it's ability to handle pca, umap, and seurat clustering.
#'
#' @param idx The names of the datasets
#' @param seed The random seed that you would like to use on the dataset
#' @return a result list with the dataset ID and the cluster it belongs to
#' @export
gficf_workflow_seurat <- function(idx, seed = 1)
{
    set.seed(seed)
    dataList <- extract_datasets(idx)
    #dataList <- extract_common_genes(dataList)
    dataList <- run_gficf(dataList)
    #dataList <- preprocess(dataList)
    dataList <- annotate_datasets(dataList)
    dataList <- select_hvg(dataList)
    dataList <- run_pca(dataList)
    dataList <- run_umap(dataList)
    dataList <- run_cluster(dataList)
    result <- cbind(dataList[[1]]@meta.data$ID[1],
                    length(dataList[[1]]@meta.data$seurat_clusters),
                    length(unique(dataList[[1]]@meta.data$seurat_clusters)))
    return(result)
}


#' Log Workflow Kmeans
#'
#' This workflow is designed to apply the log data transformation on the datasets
#' and test it's ability to handle pca, umap, and kmeans clustering.
#'
#' @param idx The names of the datasets
#' @param seed The random seed that you would like to use on the dataset
#' @return a result list with the dataset ID and the cluster it belongs to
#' @export
log_workflow_kmeans <- function(idx, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(idx)
  print("Extracted")
  #dataList <- extract_common_genes(dataList)
  dataList <- preprocess(dataList)
  print("Preprocessed")
  dataList <- annotate_datasets(dataList)
  print("Annotated")
  dataList <- run_log(dataList)
  print("Log Transformed")
  dataList <- select_hvg(dataList)
  print("HVG Selected")
  dataList <- scale_data(dataList)
  print("Data Scaled")
  dataList <- run_pca(dataList)
  print("PCA Completed")
  dataList <- run_umap(dataList)
  print("UMAP Completed")
  dataList <- run_kmeans(dataList, reduction_choosen = "pca")
  print("KMeans PCA Run")
  dataList <- run_kmeans(dataList, reduction_choosen = "umap")
  print("KMeans UMAP run")
  print(length(unique(dataList[[1]]@meta.data$kmeans_clusters_pca)))
  result <- cbind(dataList[[1]]@meta.data$ID[1],
                  length(dataList[[1]]@meta.data$kmeans_clusters_pca),
                  length(unique(dataList[[1]]@meta.data$kmeans_clusters_pca)),
                  length(unique(dataList[[1]]@meta.data$kmeans_clusters_umap)))
  return(result)
}
#' Test Log Workflow
#'
#' This function is designed to apply the log data transformation workflow on the datasets
#' and test it's ability to handle pca, umap, and clustering.
#'
#' @param idx The names of the datasets
#' @export
test_log_workflow <- function()
{
  results = cbind('idx', 'ncells', 'nclust_kmeans_pca_log', 'nclust_kmeans_umap_log')
  for (idx in datasets)
  {
    results <- rbind(results, log_workflow_kmeans(idx))
  }
  return(results)
}
#' Test GFICF Workflow
#'
#' This function is designed to apply the GFICF data transformation workflow on the datasets
#' and test it's ability to handle pca, umap, and clustering.
#'
#' @param idx The names of the datasets
#' @export
test_gficf_workflow <- function()
{
  results = cbind('idx', 'ncells', 'nclust_gficf')
  for (idx in datasets)
  {
    results <- rbind(results, gficf_workflow_seurat(idx))
  }
  return(results)
}





#' Build Seurat Pipeline with the Columns Permuted
#'
#' This functions builds the Seurat Pipeline, with experimental condition. The initial data set is appended with a the same data list,
#' with permuted columns. This allows the testing of the algorithmic stability of the Seurat Pipeline for the specified data list. The
#' Normalized information Distance is returned. If there is true algorithmic stability this score would be 1
#'
#' @param idx A data set with genes as rows and cells as columns
#' @return Adjusted Rand Index
build_seurat_columns <- function(idx, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(idx)
  dataList2 <- extract_datasets(idx)
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
#' @param idx A data set with genes as rows and cells as columns
#' @return Adjusted Rand Index
build_seurat_rows <- function(idx, seed = 1)
{
  set.seed(seed)
  dataList <- extract_datasets(idx)
  dataList2 <- extract_datasets(idx)
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
#' the rows of a new matrix, with the idx and the ARI recorded
#' @return a results matrix with ARI and idx number
#' @export
run_seurat_columns <- function()
{
  results = cbind('idx', 'ARI')
  for (idx in datasets)
  {
    result <- build_seurat_columns(idx)
    results <- rbind(results, cbind(idx, result))
  }
  return(results)
}

#' Run the Seurat Pipeline with row Permutations
#'
#' This function runs the built Seurat pipeline with row permutations on a data list. The results are appended into
#' the rows of a new matrix, with the idx and the ARI recorded
#' @return a results matrix with ARI and idx number
#' @export
run_seurat_rows <- function()
{
  results = cbind('idx', 'ARI')
  for (idx in datasets)
  {
    result <- build_seurat_rows(idx)
    results <- rbind(results, cbind(idx, result))
  }
  return(results)
}


#' Run Integration  (Need a better name)
#'
#' This functions completes an integration pipeline, using
#' 4 integration techniques. It implements the Harmony, FastMNN, CCA, and scTranform integration techniques.
#' These are implemented via the Seurat and Harmony packages. See ?run_harmony, ?run_cca, ?run_fastmnn, and ?run_sctransform
#' for more information on these pipeline structures. Each of these integration will write the integrated data set to a file
#' in the current directory
#'
#' @param idx A list of data set names that you would like to integrate with
#' the 4 pipelines
#' @export
run_integration <- function(idx)
{
  # Harmony
  data <- run_harmony(idx, batch_name = "SID")
  write_output(data[1], 'harmony')

 #FastMnn
  data <- run_fastmnn(idx, batch_column)
  write_output(data, 'fastmnn')

  # CCA
  data <- run_cca(idx)
  write_output(data, 'cca')

  # ScTransform
  data <- run_sctransform(idx)
  write_output(data, 'sctransform')
}

#Work Out Bug In Duplicate Dataset
#' Run Integration Duplicated (Need a better name)
#'
#' This functions completes an integration pipeline on a dupilcated set of the original data, using
#' 4 integration techniques. It implements the Harmony, FastMNN, CCA, and scTranform integration techniques.
#' These are implemented via the Seurat and Harmony packages. See ?run_harmony, ?run_cca, ?run_fastmnn, and ?run_sctransform
#' for more information on these pipeline structures. Each of these integration will write the integrated data set to a file
#' in the current directory
#'
#' @param idx A list of data set names that you would like to integrate with
#' the 4 pipelines
#' @param dup The number of duplications that would like to execute
#' @export
run_duplicate_integrations <- function(idx, dup)
{
  # Harmony
  dataList <- duplicate_datasets(idx, dup)
  data <- run_harmony(data, batch_column)
  write_output(data, 'harmony')

  # FastMNN
  dataList <- duplicate_datasets(idx, dup)
  data <- run_fastmnn(data, batch_column)
  write_output(data, 'fastmnn')

  # CCA
  dataList <- duplicate_datasets(idx, dup)
  data <- run_cca(dataList)
  write_output(data, 'cca')

  # ScTransform
  dataList <- duplicate_datasets(idx, dup)
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




