# scrnabench
A benchmarking tools to test downstream analysis on benchmark datasets. 'scrnabench' aims to provide all the tools necessary to test down stream analysis methods for scRNA-seq datasets, 
either on a set of benchmark data, or on the users own datasets. scrnabench makes use of the Seurat Object and several Seurat methods as well as various other downstream methods for convience of benchmarking scRNA-seq data.  


## Installing From Github

In the R terminal execute the following command 

```
devtools::install_github("NWhitener/scrnabench")
```
This will install the package. To use the pasckage, call the package like any other R package 

```
library(scrnabench)
```

## Common Issues 

Below are solutions to some of the common issues while using the package. 

### Vector Memory Exhausted 

This issue may occur when using the workflows if R does not have access to enough memory. To fix this execute the following commands. 

```
library(usethis) 
usethis::edit_r_environ()
```

This will open the R environment file. Add the following command, with your request memory size and restart your R session to increase the memory and 
complete the workflows 

```
R_MAX_VSIZE=100Gb
```

# Quick Start Guide 

## Dataset Download
In order to download the dataset use the download_data function. This will download the dataset from [Zenodo](https://zenodo.org/record/6617997).
To download the data provide a path for the download. This function will download and automatically load the gene_counts.RDS file which can be used for benchmarking. 
This function only needs to be run once. 

```
datasets = dowload_data(path = "/Users/Downloads")  
```

## Dataset Load 
If you have already downloaded the data or would like to use the demo dataset use the load_data() function. This will load the dataset into memory for usage. To 
load the full dataset, provide the path of the gene_counts.RDS file from the download.

```
#Load the demo dataset 
datasets = load_data(demo = TRUE) 

#Load the full dataset 
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
```

## Workflows 
Now that the data is loaded, the clustering and integration workflows are ready to use. Below we show the syntax for the workflows.

### Clustering

The clustering workflow completes the full clustering pipeline. The log transformation includes preprocessing, transformation, Highly Variable Gene Selection, 
scaling, dimensionality reduction, and clustering. The tfidf transformation includes transformation, Highly Variable Gene Selection, scaling, dimensionality reduction and 
clustering. The cluster method can be KMeans Clustering or Seurat's Phenograph clustering. 

```
#Clustering Workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_clustering_workflow(dataList, method = 'kmeans', transformationType = 'log', seed = 17, numberClusters = 10)
```

### Harmony Integration 

The harmony integration workflow completes the full harmony pipeline. This includes extracting common genes, merging datasets, preprocessing, transformation,
Highly Variable Gene selection, scaling, integration via harmony, dimensionality reduction, and clustering. The cluster method can be KMeans Clustering or Seurat's Phenograph clustering.

```
#Harmony workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_harmony_integration_workflow(dataList, method = 'kmeans', seed = 17, numberClusters = 10)
```

### FastMnn Integration 

The fastmnn integration workflow completes the full fastmnn pipeline. This includes extracting common genes, merging datasets, preprocessing, transformation,
Highly Variable Gene selection, integration via fastmnn, dimensionality reduction, and clustering. The cluster method can be KMeans Clustering or Seurat's Phenograph clustering.

```
#FastMNN workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_fastmnn_integration_workflow(dataList, method = 'kmeans', seed = 17, numberClusters = 10)
```
### CCA Integration 

The cca integration workflow completes the full cca pipeline. This includes preprocessing, transformation,
Highly Variable Gene selection, integration via cca, scaling, dimensionality reduction, and clustering. The cluster method can be KMeans Clustering or Seurat's Phenograph clustering.

```
#cca workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_cca_integration_workflow(dataList, method = 'kmeans', seed = 17, numberClusters = 10)
```

### Sctransform Integration 

The SCtransform integration workflow completes the full cca pipeline. This includes preprocessing, integration via sctransform scaling, dimensionality reduction, and clustering. 
The cluster method can be KMeans Clustering or Seurat's Phenograph clustering.

```
#SCtransform workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_sctransform_integration_workflow(dataList, method = 'kmeans', seed = 17, numberClusters = 10)
```
