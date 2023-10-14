# scrnabench
'scrnabench' provides tools for the creation of reference and metamorphic datasets for benchmarking of scRNA-seq data analysis methods.
This R package focuses on metamorphic stability of cluster analysis and data integration tools. 

## Installing From Github

The `BiocManager` is needed for installing scrnabench. Execute the following commands:
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("batchelor")
```

If you do not have `devtools` installed, run this command:

```
install.packages("devtools")
```

In the R terminal, execute the following command to install the package:
```
devtools::install_github("NWhitener/scrnabench")
```

To use the package:
```
library(scrnabench)
```

To download the raw data, use the **download_data()** function. This will download 48 raw reference datasets from [Zenodo](https://zenodo.org/record/6617997). The datasets are downloaded as one R object called 'gene_counts.RData'.
To save the dataset, provide a filepath to the desired location, for example, **path = "/User/Downloads"**.  This function will download and automatically load the raw datasets into scrnabench. This function only needs to be run once.                                                                          

```
datasets = download_data(path = "/Users/Downloads")
```

## Common Issues 

Below are some of the solutions to common issues that may arise when installing or using the package. 

### ERROR: lazy loading failed for package ‘SeuratWrappers’

This issue may occur during the installation. To fix this issue, execute the following commands:

```
remotes::install_github('satijalab/seurat-wrappers')
devtools::install_github("NWhitener/scrnabench")
```

### API rate limit exceeded

This issue may occur when you call a GitHub API more than 5,000 times within a 60-minute window, even if you are authenticated. To fix this issue, see: [post](https://gist.github.com/Z3tt/3dab3535007acf108391649766409421). 

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

## Dataset Load 
If you have already downloaded the raw data to a directory on your local computer or would like to use the demo dataset, call the **load_data()** function.  Provide the path of the downloaded data file.  This will load the dataset for benchmarking. This function should be used anytime the data needs to be loaded.

```
#Load the demo dataset 
datasets = load_data(demo = TRUE, path = "/Users/Downloads") 

#Load the full dataset 
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
```

The demo file is smaller than the full dataset, comprising of 2 raw reference datasets, and can be used for rapid testing and experimentation.

## Workflows 
Below we have provided sample workflows to demonstrate the functionality of the scrnabench package.

### Cluster Analysis 

The clustering workflow performs the baseline cluster analysis pipeline. The following steps are performed: filtering, annotation, data transformation, Highly Variable Gene (HVG) selection, 
scaling, dimensionality reduction, and clustering. The package supports 2 baseline clustering methods, Kmeans clustering and [Seurat's](https://satijalab.org/seurat/) Phenograph clustering. 

For example, to cluster the demo datasets into 10 clusters using Kmeans algorithm, the following command can be used: 

```
#Clustering Workflow 
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_clustering_workflow(dataList, method = 'kmeans', transformationType = 'log', seed = 1, numberClusters = 10)
```
To cluster the demo datasets using Seurat's Phenograph, do the following: 
 ```
 #Clustering Workflow 
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_clustering_workflow(dataList, method = 'seurat', transformationType = 'log', seed = 1)
 ```

### Metamorphic Benchmarking

The metamorphic benchmarking workflow performs metamorphic cluster analysis. This workflow perturbs raw datasets using up to 6 metamorphic relations, completes preprocessing and cluster analysis. See [protocols](https://core.ac.uk/download/pdf/228032568.pdf) for details about metamorphic testing. The benchmarking report summarizes clustering changes of reference and metamorphic datasets. 

For example, to cluster the data into 10 clusters using Kmeans algorithm, the following command can be used: 

```
#Metamorphic Testing Workflow 
datasets = load_data(demo = F, path = 'Users/Downloads') 
dataList = extract_datasets(datasets) 
MetamorphicReport = run_metamorphic_test_workflow(dataList, method = 'kmeans', transformation = "log", metamorphicTests= c(6,5,4,3,2,1))
```
To cluster the data using Seurat's Phenograph, do the following: 
 ```
datasets = load_data(demo = F, path = '/Users/Downloads') 
dataList = extract_datasets(datasets) 
MetamorphicReport = run_metamorphic_test_workflow(dataList, method = 'seurat', transformation = "log", metamorphicTests= c(6,5,4,3,2,1))
```

### Harmony Integration 

The Harmony integration workflow performs data integration using the [harmony](https://github.com/immunogenomics/harmony) package. The workflow steps include: extraction of common genes, dataset merging, filtering, annotation, data transformation,
HVG selection, scaling, harmonization, dimensionality reduction, and clustering. Two baseline clustering algorithms are supported, Kmeans and Seurat Phenograph.


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


### fastMNN Integration 

The fastMNN integration workflow completes the data integration pipeline, using the SeuratWrapper's [fastMNN](https://github.com/satijalab/seurat-wrappers) functions (LINK). Here the steps include: extraction of common genes, dataset merging, filtering, annotation, data transformation,
HVG selection, integration via fastmnn, dimensionality reduction, and clustering. The cluster method can be KMeans Clustering or Seurat's Phenograph clustering.

To integrate the data and complete Kmeans clustering, the following command can be used: 
```
#FastMNN workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_fastmnn_integration_workflow(dataList, method = 'kmeans', seed = 1, numberClusters = 10)
```
To integrate the data and complete Seurat's Phenograph, do the following: 
```
#FastMNN workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_fastmnn_integration_workflow(dataList, method = 'seurat', seed = 1)
```


### CCA Integration 

The cca integration workflow completes the data integration usign Seurat's CCA pipeline.  Each dataset is preproceesed separately, then integrated via CCA, followed by scaling, dimensionality reduction, and clustering. The preprocess consists of filtering, annotation, data transformation, and HVG selection. Two baseline clustering algorithms are supported, Kmeans and Seurat Phenograph.


To integrate the data and complete Kmeans clustering, the following command can be used: 
```
#cca workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_cca_integration_workflow(dataList, method = 'kmeans', seed = 1, numberClusters = 10)
```
To integrate the data and complete Seurat's Phenograph, do the following: 
```
#cca workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_cca_integration_workflow(dataList, method = 'seurat', seed = 1)
```

### sctransform Integration 

The sctransform integration workflow completes the Seurat scransform pipeline. This includes filtering, annotation, integration via sctransform, scaling, dimensionality reduction, and clustering. Two baseline clustering algorithms are supported, Kmeans and Seurat Phenograph.


To integrate the data and complete Kmeans clustering, the following command can be used: 
```
#SCtransform workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_sctransform_integration_workflow(dataList, method = 'kmeans', seed = 1, numberClusters = 10)
```
To integrate the data and complete Seurat's Phenograph, do the following: 
```
#SCtransform workflow
datasets = load_data(demo = FALSE, path = "/Users/Downloads")
dataList = extract_datasets(datasets)
dataList = run_sctransform_integration_workflow(dataList, method = 'Seurat', seed = 1)
```

### Publications
The srnabench has been used in the following publications.

```
@inproceedings{zhao2023ensemble,
  title={An Ensemble Machine Learning Approach for Benchmarking and Selection of scRNA-seq Integration Methods},
  author={Zhao, Konghao and Bhandari, Sapan and Whitener, Nathan P and Grayson, Jason M and Khuri, Natalia},
  booktitle={Proceedings of the 14th ACM International Conference on Bioinformatics, Computational Biology, and Health Informatics},
  pages={1--10},
  year={2023}
}

@Article{jpm13020183,
AUTHOR = {Zhao, Konghao and Grayson, Jason M. and Khuri, Natalia},
TITLE = {Multi-Objective Genetic Algorithm for Cluster Analysis of Single-Cell Transcriptomes},
JOURNAL = {Journal of Personalized Medicine},
VOLUME = {13},
YEAR = {2023},
NUMBER = {2},
ARTICLE-NUMBER = {183},
URL = {https://www.mdpi.com/2075-4426/13/2/183},
PubMedID = {36836417},
ISSN = {2075-4426},
DOI = {10.3390/jpm13020183}
}
``` 
