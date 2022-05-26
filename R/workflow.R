run_cca <- function(data.list)
{
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })

  features <- Suerat::SelectIntegrationFeatures(object.list = data.list)
  k.filter <- min(200, min(sapply(data.list, ncol)))
  data.anchors <- Seurat::FindIntegrationAnchors(object.list = data.list, anchor.features = features, k.filter=k.filter)
  data.combined <- Seurat::IntegrateData(anchorset = data.anchors)
  DefaultAssay(data.combined) <- "integrated"

  data.combined <- Seurat::ScaleData(data.combined, verbose = FALSE)
  data.combined <- Seurat::RunPCA(data.combined, npcs = 30, verbose = FALSE)
  data.combined <- Seurat::RunUMAP(data.combined, reduction = "pca", dims = 1:30)
  data.combined <- Seurat::FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
  data.combined <- Seurat::FindClusters(data.combined, resolution = 0.5)

  return(data.combined)
}




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


run_fastmnn <- function(data, batch_name)
{
  data <- Seurat::FindVariableFeatures(data, nfeatures=2000)
  data <- SeuratWrapper::RunFastMNN(object.list = SplitObject(data, split.by = batch_name))
  data <- Seurat::RunUMAP(data, reduction = "mnn", dims = 1:30) %>%
    Seurat::FindNeighbors(reduction = "mnn", dims = 1:30) %>%
    Seurat::FindClusters()
  DefaultAssay(data) <- "mnn.reconstructed"
  return(data)

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


#' Run Seurat Pipeline with the Columns Permuted
#'
#' This functions runs the Seurat Pipeline, with experimental condition. The initial data list is appended with a the same data list,
#' with permuted columns. This allows the testing of the algorithmic stability of the Seurat Pipeline for the specified data list. The
#' Normalized information Distance is returned. If there is true algorithmic stability this score would be 1
#'
#' @param idx A data list with genes as rows and cells as columns
#' @return Normalized Information Distance
#' @export
#'
run_seurat_columns <- function(idx)
{
  data.list <- extract_datasets(idx)
  data.list2 <- extract_datasets(idx)
  data.list2 <- permute_columns(data.list2)
  data.list <- append(data.list, data.list2)
  data.list <- preprocess_data(data.list)
  data.list <- run_log(data.list)
  data.list <- run_seurat(data.list)

  x = merge(data.list[[1]]@meta.data, data.list[[2]]@meta.data, by="row.names", all=TRUE)
  return(aricode::NID(x$seurat_clusters.x, x$seurat_clusters.y))
}

#' Run Seurat Pipeline with the Rows Permuted
#'
#' This functions runs the Seurat Pipeline, with experimental condition. The initial data list is appended with a the same data list,
#' with permuted rows. This allows the testing of the algorithmic stability of the Seurat Pipeline for the specified data list. The
#' Normalized information Distance is returned. If there is true algorithmic stability this score would be 1
#'
#' @param idx A data list with genes as rows and cells as columns
#' @return Adjusted Rand Index
#' @export
#'
build_seurat_rows <- function(idx, seed = 1 )
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






run_sctransform <- function(data.list)
{
  data.list <- lapply(X = data.list, FUN = SCTransform)
  features <- Seurat::SelectIntegrationFeatures(object.list = data.list, nfeatures=2000)
  data.list <- Seurart::PrepSCTIntegration(object.list = data.list, anchor.features = features)
  k.filter <- min(200, min(sapply(data.list, ncol)))
  data.anchors <- Seurat::FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT", anchor.features = features, k.filter = k.filter)
  data.combined <- Seurat::IntegrateData(anchorset = data.anchors, normalization.method = "SCT")
  DefaultAssay(data.combined) <- "integrated"

  data.combined <- Seurat::ScaleData(data.combined, verbose = FALSE)
  data.combined <- Seurat::RunPCA(data.combined, npcs = 30, verbose = FALSE)
  data.combined <- Seurat::RunUMAP(data.combined, reduction = "pca", dims = 1:30)
  data.combined <- Seurat::FindNeighbors(data.combined, reduction = "pca", dims = 1:30)
  data.combined <- Seurat::FindClusters(data.combined, resolution = 0.5)

  return(data.combined)

}


run_workflow <- function(idx)
{
  # Harmony
  data <- run_harmony(data, batch_column)
  write_output(data[1], 'harmony')

  # FastMNN
  data.list <- extract_datasets(idx)
  data.list <- extract_common_genes(data.list)
  data.list <- merge_datasets(data.list, intersect=TRUE)
  #data.list <- run_gficf(data.list)
  data.list <- preprocess_data(data.list)
  data.list <- run_log(data.list) #LOG
  #data <- annotate_seurat_object(data.list[[1]])
  data <- run_fastmnn(data, batch_column)
  write_output(data, 'fastmnn')

  # CCA
  data.list <- extract_datasets(idx)
  data.list <- preprocess_data(data.list)
  #data.list <- run_gficf(data.list)
  data.list <- run_log(data.list)
  data <- run_cca(data.list)
  #data <- annotate_seurat_object(data)
  write_output(data, 'cca')

  # ScTransform
  data.list <- extract_datasets(idx)
  data.list <- preprocess_data(data.list)
  data <- run_sctransform(data.list)
  #data <- annotate_seurat_object(data)
  write_output(data, 'sctransform')
}


run_workflow2 <- function(idx, dup)
{
  # Harmony
  data.list <- duplicate_datasets(idx, dup)
  data.list <- extract_common_genes(data.list)
  data.list <- merge_datasets(data.list, intersect=TRUE)
  #data.list <- run_gficf(data.list)
  data.list <- preprocess_data(data.list)
  data.list <- run_log(data.list) #LOG

 # data <- annotate_seurat_object(data.list[[1]])
  data <- run_harmony(data, batch_column)
  write_output(data, 'harmony')

  # FastMNN
  data.list <- duplicate_datasets(idx, dup)
  data.list <- extract_common_genes(data.list)
  data.list <- merge_datasets(data.list, intersect=TRUE)
  #data.list <- run_gficf(data.list)
  data.list <- preprocess_data(data.list)
  data.list <- run_log(data.list) #LOG
 # data <- annotate_seurat_object(data.list[[1]])
  data <- run_fastmnn(data, batch_column)
  write_output(data, 'fastmnn')

  # CCA
  data.list <- duplicate_datasets(idx, dup)
  data.list <- preprocess_data(data.list)
  #data.list <- run_gficf(data.list)
  data.list <- run_log(data.list)
  data <- run_cca(data.list)
  #data <- annotate_seurat_object(data)
  write_output(data, 'cca')

  # ScTransform
  data.list <- duplicate_datasets(idx, dup)
  data.list <- preprocess_data(data.list)
  data <- run_sctransform(data.list)
  #data <- annotate_seurat_object(data)
  write_output(data, 'sctransform')
}


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


