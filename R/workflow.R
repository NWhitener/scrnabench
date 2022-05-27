#' Run CCA
#'
#' This function runs the data integration protocol detailed in the Seurat "Introduction to scRNA-seq
#' integration" found at https://satijalab.org/seurat/articles/integration_introduction.html. This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param data.list A data list of data sets to integrate using the cca protocol
#' @return data.combined A data list of the combined data from the cca protocol
#' @export
run_cca <- function(data.list)
{

  data.list <- extract_datasets(idx)
  data.list <- preprocess_data(data.list)
  #data.list <- run_gficf(data.list)
  data.list <- run_log(data.list)

  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- select_hvg(x)
  })

  features <- Suerat::SelectIntegrationFeatures(object.list = data.list)
  k.filter <- min(200, min(sapply(data.list, ncol)))
  data.anchors <- Seurat::FindIntegrationAnchors(object.list = data.list, anchor.features = features, k.filter=k.filter)
  data.combined <- Seurat::IntegrateData(anchorset = data.anchors)
  DefaultAssay(data.combined) <- "integrated"

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
#' @param data.list A data list of data sets to integrate using the harmony protocol
#' @return data.combined A data list of the combined data from the harmony protocol
#' @export
run_harmony <- function(idx, batch_name)
{
    data.list <- extract_datasets(idx)
    data.list <- extract_common_genes(data.list)
    data.list <- merge_datasets(data.list, intersect=TRUE)
    #data.list <- run_gficf(data.list)
    data.list <- preprocess_data(data.list)
    data.list <- run_log(data.list) #LOG
    data.list <- select_hvg(data.list)
    data.list = scale_data(data.list)
    data.list = run_pca(data.list)
    data.list = harmony::RunHarmony(group.by.vars = batch_name)
    data.list = run_umap(data.list, "harmony")
    data.list = run_cluster(data.list, "harmony")
    retunr(data.list)
}

#' Run fastmnn
#'
#' This function runs the data integration protocol detailed in the Seurat ## FIND VING. AND LINK . This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param data.list A data list of data sets to integrate using the fastmnn protocol
#' @return data.combined A data list of the combined data from the fastmnn protocol
#' @export
run_fastmnn <- function(idx, batch_name)
{

  data.list <- extract_datasets(idx)
  data.list <- extract_common_genes(data.list)
  data.list <- merge_datasets(data.list, intersect=TRUE)
  data.list <- preprocess_data(data.list)
  data.list <- run_log(data.list) #LOG
  #data.list<- run_gficf(data.list)
  data.list <- select_hvg(data.list)
  data.list <- SeuratWrapper::RunFastMNN(object.list = SplitObject(data.list, split.by = batch_name))
  data.list <- run_umap(data.list)
  data.list <- run_cluster(data.list)
  DefaultAssay(data.list) <- "mnn.reconstructed"
  return(data.list)

}


#' Run sctransform
#'
#' This function runs the data integration protocol detailed in the Seurat "Using sctransform in Seurtat"
#' found at https://satijalab.org/seurat/articles/sctransform_vignette.html . This function completes the entire
#' protocol on either gficf or log transformed data.
#'
#' @param data.list A data list of data sets to integrate using the sctransform protocol
#' @return data.combined A data list of the combined data from the sctransform protocol
#' @export
run_sctransform <- function(data.list)
{

  data.list <- extract_datasets(idx)
  data.list <- preprocess_data(data.list)
  data.list <- lapply(X = data.list, FUN = SCTransform)
  features <- Seurat::SelectIntegrationFeatures(object.list = data.list, nfeatures=2000)
  data.list <- Seurart::PrepSCTIntegration(object.list = data.list, anchor.features = features)
  k.filter <- min(200, min(sapply(data.list, ncol)))
  data.anchors <- Seurat::FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = features, k.filter = k.filter)
  data.combined <- Seurat::IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
  DefaultAssay(data.combined) <- "integrated"
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
#' @param data.list A data list of Data sets in the Gene in Row and Cell in Columns format
#' @export
run_seurat <- function(data.list)
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- RunPCA(x, npcs = 30, verbose = FALSE)
    x <- RunUMAP(x, reduction = "pca", dims = 1:30)
    x <- FindNeighbors(x, reduction = "pca", dims = 1:30)
    x <- FindClusters(x, resolution = 0.5)
  })
  return(data.list)
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
  data.list <- extract_datasets(idx)
  data.list2 <- extract_datasets(idx)
  data.list2 <- permute_columns(data.list2)
  data.list <- append(data.list, data.list2)
  data.list <- preprocess_data(data.list)
  data.list <- run_log(data.list)
  data.list <- run_seurat(data.list)
  x = merge(data.list[[1]]@meta.data, data.list[[2]]@meta.data, by="row.names", all=TRUE)
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
  data.list <- extract_datasets(idx)
  data.list2 <- extract_datasets(idx)
  data.list2 <- permute_rows(data.list2)
  data.list <- append(data.list, data.list2)
  data.list <- preprocess_data(data.list)
  data.list <- run_log(data.list)
  data.list <- run_seurat(data.list)
  x = merge(data.list[[1]]@meta.data, data.list[[2]]@meta.data, by="row.names", all=TRUE)
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
  data <- run_harmony(data.list, batch_column)
  write_output(data[1], 'harmony')

 #FastMnn
  data <- run_fastmnn(data.list, batch_column)
  write_output(data, 'fastmnn')

  # CCA
  data <- run_cca(data.list)
  write_output(data, 'cca')

  # ScTransform
  data <- run_sctransform(data.list)
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
  data.list <- duplicate_datasets(idx, dup)
  data <- run_harmony(data, batch_column)
  write_output(data, 'harmony')

  # FastMNN
  data.list <- duplicate_datasets(idx, dup)
  data <- run_fastmnn(data, batch_column)
  write_output(data, 'fastmnn')

  # CCA
  data.list <- duplicate_datasets(idx, dup)
  data <- run_cca(data.list)
  write_output(data, 'cca')

  # ScTransform
  data.list <- duplicate_datasets(idx, dup)
  data <- run_sctransform(data.list)

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
  print(dim(as.data.frame(GetAssayData(data))[VariableFeatures(data),]))
  path = paste("./", prefix, sep="")

  write.csv(t(as.data.frame(GetAssayData(data))[VariableFeatures(data),]), paste(path, "_hvg.csv", sep=""), quote=F)
  write.csv(Embeddings(data, 'umap'), paste(path, "_umap.csv", sep=""), quote=F)
  if(prefix=='harmony')
  {
    write.csv(Embeddings(data, 'harmony'), paste(path, "_harmonized.csv", sep=""), quote=F)
  }
  else if (prefix == 'fastmnn')
  {
    write.csv(Embeddings(data, 'mnn'), paste(path, "_harmonized.csv", sep=""), quote=F)
  }
  else
  {
    write.csv(Embeddings(data, 'pca'), paste(path, "_pca.csv", sep=""), quote=F)
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


