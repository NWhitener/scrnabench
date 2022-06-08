run_harmony3 <- function(data, batch_name)
{
  data <- Seurat::FindVariableFeatures(data, nfeatures = 2000) %>%
    Seurat::ScaleData() %>%
    Seurat::RunPCA(verbose = FALSE) %>%
    Seurat::RunHarmony(group.by.vars = batch_name) %>%
    Seurat::RunUMAP(reduction = "harmony", dims = 1:30) %>%
    Seurat::FindNeighbors(data, reduction = "harmony", dims = 1:30) %>%
    Seurat::FindClusters()
}



run_workflow <- function(idx)
{
  # Harmony
  dataList <- extract_datasets(idx)
  dataList <- extract_common_genes(dataList)
  dataList <- merge_datasets(dataList)
  #dataList <- run_gficf(dataList)
  dataList <- preprocess(dataList)
  data <- annotate_datasets(dataList[[1]])
  dataList <- run_log(dataList) #LO
  dataList <- run_
  data <- run_harmony3(data, batch_column)
  #write_output(data, 'harmony')
}
