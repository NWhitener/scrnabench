# scrnabench
'scrnabench' aims to provide the tools necessary to test downstream analysis methods on 48 reference scRNA-seq datasets.
This R package focuses on metamorphic stability of clustert analysis and data integration. 

## Installing From Github

In the R terminal execute the following command 

```
devtools::install_github("NWhitener/scrnabench")
```
This will install the package. To use the package:
```
library(scrnabench)
```

## Common Issues 

Below are some of the solutions to some of the common issues that may appear while using the package. 

### Vector Memory Exhausted 

This issue may occur when R does not have a large enough max vector heap size. To fix this issue, execute the following commands: 

```
library(usethis) 
usethis::edit_r_environ()
```

This will open the R environment file. Add the following command, requesting larger memory (example 100 Gb) and restart your R session. 

```
R_MAX_VSIZE=100Gb
```

# Quick Start Guide 

## Dataset Download
In order to download the dataset use the **download_data()** function. This will download the dataset from [Zenodo](https://zenodo.org/record/6617997).
To download the data provide a path to the desired directory where the data should be downloaded. For example, **path = "/User/Downloads"**.  This function will download and automatically load the data file which can be immediately used for benchmarking. This function only needs to be run once. 

```
datasets = dowload_data(path = "/Users/Downloads")  
```

## Dataset Load 
If you have already downloaded the data to a directory on your local computer or would like to use the demo dataset, call the **load_data()** function.  Provide the path of the downloaded data file.  This will load the dataset for benchmarking. This function should be used anytime the data needs to be loaded.

```
#Load the demo dataset 
datasets = load_data(demo = TRUE, path = "/Users/Downloads") 

#Load the full dataset 
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
```

The demo dataset is smaller than the full dataset, comprising of 2 data files, and can be used for rapid testing and experimentation.

## Workflows 
Below we have provided sample workflows to demonstrate the functionality of the scrnabench package

### Cluster Analysis 

The clustering workflow completes the standard clustering pipeline. The following steps are performed: filtering, annotation, data transformation, Highly Variable Gene selection(HVG), 
scaling, dimensionality reduction, and clustering. The pacakge supports two baseline clustering methods, Kmeans clustering and Seurat's Phenograph clustering. 

For example, to cluster the data into 10 clusters using Kmeans algorithm, the following command can be used: 

```
#Clustering Workflow 
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_clustering_workflow(dataList, method = 'kmeans', transformationType = 'log', seed = 1, numberClusters = 10)
```
To cluster the data using Seurat's Phenograph, do the following: 
 ```
 #Clustering Workflow 
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_clustering_workflow(dataList, method = 'seurat', transformationType = 'log', seed = 1)
 ```

### Harmony Integration 

The Harmony integration workflow completes a data integration pipeline using the harmony pacakge. The workflow steps include, extraction of common genes, dataset merging, preprocessing, data transformation,
HVG, scaling, harmonization, dimensionality reduction, and clustering. Two baseline clustering algorithms are supported, Kmeans and Seurat Phenograph


For example, to integrate the data and complete Kmeans clustering, the following command can be used: 
```
#Harmony workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_harmony_integration_workflow(dataList, method = 'kmeans', seed = 1, numberClusters = 10)
```
To integrate the data and complete Seurat's Phenograph, do the following: 

```
#Harmony workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_harmony_integration_workflow(dataList, method = 'seurat', seed = 1)
```


### FastMnn Integration 

The fastmnn integration workflow completes the full fastmnn pipeline. This includes extracting common genes, merging datasets, preprocessing, transformation,
Highly Variable Gene selection, integration via fastmnn, dimensionality reduction, and clustering. The cluster method can be KMeans Clustering or Seurat's Phenograph clustering.

```
#FastMNN workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_fastmnn_integration_workflow(dataList, method = 'kmeans', seed = 1, numberClusters = 10)
```
### CCA Integration 

The cca integration workflow completes the full cca pipeline. This includes preprocessing, transformation,
Highly Variable Gene selection, integration via cca, scaling, dimensionality reduction, and clustering. The cluster method can be KMeans Clustering or Seurat's Phenograph clustering.

```
#cca workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_cca_integration_workflow(dataList, method = 'kmeans', seed = 1, numberClusters = 10)
```

### Sctransform Integration 

The SCtransform integration workflow completes the full cca pipeline. This includes preprocessing, integration via sctransform scaling, dimensionality reduction, and clustering. 
The cluster method can be KMeans Clustering or Seurat's Phenograph clustering.

```
#SCtransform workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_sctransform_integration_workflow(dataList, method = 'kmeans', seed = 1, numberClusters = 10)
```
