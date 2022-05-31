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
  data.list <- extract_datasets(idx)
  data.list <- extract_common_genes(data.list)
  data.list <- merge_datasets(data.list, intersect=TRUE)
  #data.list <- run_gficf(data.list)
  data.list <- preprocess_data(data.list)
  data.list <- run_log(data.list) #LOG

  data <- annotate_seurat_object(data.list[[1]])
  data <- run_harmony3(data, batch_column)
  #write_output(data, 'harmony')
}
