run_harmony3 <- function(data, batch_name)
{
  data <- FindVariableFeatures(data, nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE) %>%
    RunHarmony(group.by.vars = batch_name) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(data, reduction = "harmony", dims = 1:30) %>%
    FindClusters()
}



run_workflow <- function(idx)
{
  # Harmony
  dataList <- extract_datasets(idx)
  dataList <- extract_common_genes(dataList)
  dataList <- merge_datasets(dataList, intersect=TRUE)
  #dataList <- run_gficf(dataList)
  dataList <- preprocess_data(dataList)
  dataList <- run_log(dataList) #LOG

  data <- annotate_seurat_object(dataList[[1]])
  data <- run_harmony3(data, batch_column)
  #write_output(data, 'harmony')
}
