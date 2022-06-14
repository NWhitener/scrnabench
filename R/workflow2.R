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
