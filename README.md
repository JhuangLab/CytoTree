
# CytoTree <img src="https://github.com/JhuangLab/CytoTree/blob/master/inst/figures/logo.png" align="right" height=150 width=150/>

CytoTree is an R package to implement cellular subpopulations identification, trajectory inference, pseudotime estimation and visualization for flow and mass cytometry data. This package is developed and maintained by [JhuangLab](https://github.com/JhuangLab) at Shanghai Institute of Hematology.

See the quick start tutorial of CytoTree, please visit [Quick start of CytoTree](https://ytdai.github.io/CytoTree/Quick_start.html).

See the basic tutorial of CytoTree, please visit [Tutorial of CytoTree](https://ytdai.github.io/CytoTree/basic.html).

See time-course data analysis of CytoTree, please visit [Time-course workflow of CytoTree](https://ytdai.github.io/CytoTree/Time_course.html).


Use cases could be found at: 

https://github.com/JhuangLab/CytoTree-dataset


You can view and clone the use cases of CytoTree on GitHub at by `git clone https://github.com/JhuangLab/CytoTree-dataset`


## 1 Introduction

Multidimensional single-cell-based flow and mass cytometry  enable ones to analyze multiple single-cell parameters and identify cellular populations. 
Based on classical software for analyzing [Flow Cytometry Standard](https://en.wikipedia.org/wiki/Flow_Cytometry_Standard) (FCS) data such as [`flowSOM`](https://bioconductor.org/packages/release/bioc/html/FlowSOM.html)[1] and [`SPADE`](https://github.com/nolanlab/spade)[2], methods for inferencing cellular trajectory during a biological process are very important. 
To objectively inference differential trajectory based on time courses FCS data, we present [`CytoTree`](https://github.com/JhuangLab/CytoTree), a trajectory inference and visualization toolkit for flow and mass cytometry data. 

`CytoTree` can help you to perform four main types of analysis:

- **Clustering**. `CytoTree` can help you to discover and identify subtypes of cells. 

- **Dimensionality Reduction**. Several dimensionality reduction methods are provided in `CytoTree` package such as Principal Components Analysis (PCA), t-distributed Stochastic Neighbor Embedding (tSNE), Diffusion Maps and Uniform Manifold Approximation and Projection (UMAP). CytoTree provides both cell-based and cluster-based dimensionality reduction.

- **Trajectory Inference**. `CytoTree` can help you to construct the cellular differential based on minimum spanning tree (MST) algorithm. 

- **Pseudotime and Intermediate states definition**. The root cells need to be defined by users. The trajctroy value will be calculated based on Shortest Path from root cells and leaf cells using R `igraph` package. Subset FCS data set in `CytoTree` and find the key intermediate cell states based on trajectory value.

## 2 Installation

### 2.1 From Github

This requires the `devtools` package to be installed first.

```

# If not already installed
install.packages("devtools") 
devtools::install_github("JhuangLab/CytoTree")

library(CytoTree)

```


## 3 Quick start (Standard Workflow)

``` {r}

# Loading packages
suppressMessages({
library(ggplot2)
library(CytoTree)
library(flowCore)
library(stringr)
})

# Read fcs files
fcs.path <- system.file("extdata", package = "CytoTree")
fcs.files <- list.files(fcs.path, pattern = '.FCS$', full = TRUE)

fcs.data <- runExprsMerge(fcs.files, comp = FALSE, transformMethod = "none")

# Refine colnames of fcs data
recol <- c(`FITC-A<CD43>` = "CD43", `APC-A<CD34>` = "CD34", 
           `BV421-A<CD90>` = "CD90", `BV510-A<CD45RA>` = "CD45RA", 
           `BV605-A<CD31>` = "CD31", `BV650-A<CD49f>` = "CD49f",
           `BV 735-A<CD73>` = "CD73", `BV786-A<CD45>` = "CD45", 
           `PE-A<FLK1>` = "FLK1", `PE-Cy7-A<CD38>` = "CD38")
colnames(fcs.data)[match(names(recol), colnames(fcs.data))] = recol
fcs.data <- fcs.data[, recol]

day.list <- c("D0", "D2", "D4", "D6", "D8", "D10")
meta.data <- data.frame(cell = rownames(fcs.data),
                        stage = str_replace(rownames(fcs.data), regex(".FCS.+"), "") )
meta.data$stage <- factor(as.character(meta.data$stage), levels = day.list)

markers <- c("CD43","CD34","CD90","CD45RA","CD31","CD49f","CD73","CD45","FLK1","CD38")

# Build the CYT object
cyt <- createCYT(raw.data = fcs.data, markers = markers,
                   meta.data = meta.data,
                   normalization.method = "log",
                   verbose = TRUE)

# See information
cyt

# Standard workflow of CytoTree
cyt <- runCluster(cyt)
cyt <- processingCluster(cyt)
cyt <- runFastPCA(cyt)
cyt <- runTSNE(cyt)
cyt <- runDiffusionMap(cyt)
cyt <- runUMAP(cyt)
cyt <- buildTree(cyt, dim.type = "umap", dim.use = 1:2)
cyt <- defRootCells(cyt, root.cells = 1)
cyt <- runPseudotime(cyt)
cyt <- defLeafCells(cyt, leaf.cells = 2)
cyt <- runWalk(cyt)


```

## 4 Reported bugs and solutions

If there is any error in installing or librarying the `CytoTree` package, please contact us via e-mail forlynna@sjtu.edu.cn

## 5 Version History

Jun 02, 2020
 - Version 0.99.3
 - Changes:
   - Remove if (FALSE) in examples

May 10, 2020
 - Version 0.99.0
 - Changes:
   - First commit of CytoTree

## 6 Note


The previous version of `CytoTree` is `flowSpy` **[link to GitHub](https://github.com/JhuangLab/CytoTree) and [link to Bioconductor](https://bioconductor.org/packages/flowSpy/)**. To improve the identification and avoid awkward duplication of names in some situations, we changed the name of `flowSpy` to `CytoTree`. `CytoTree` more fits the functional orientation of this software.

We apologized for the inconvenience.


## 7 Reference

[1] Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2019). FlowSOM: Using
  self-organizing maps for visualization and interpretation of cytometry data.
  http://www.r-project.org, http://dambi.ugent.be.

[2] Qiu, P., et al., Extracting a cellular hierarchy from high-dimensional cytometry data with SPADE. Nat Biotechnol, 2011. 29(10): p.886-91.





