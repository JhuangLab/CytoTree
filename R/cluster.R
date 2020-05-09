#'
#' Specific Clustering Method Toolkits
#'
#' @name runCluster
#'
#' @description Compute a specific clustering using the combined flow
#'    cytometry data. "som" \code{\link[FlowSOM]{SOM}}, "hclust" \code{\link[stats]{hclust}},
#'    "clara" \code{\link[cluster]{clara}}, "phenograph", "kmeans" \code{\link[stats]{kmeans}} are
#'    provided.
#'
#' @param object an FSPY object
#' @param cluster.method character. Four clustering method are provided: som, clara, kmeans and phenograph.
#'    Clustering method "hclust" and "mclust" are not recommended because of long computing time.
#' @param verbose logic. Whether to print calculation progress.
#' @param ... options to pass on to the clustering functions.
#'
#' @seealso \code{\link[FlowSOM]{SOM}}, \code{\link[stats]{hclust}},
#'    \code{\link[cluster]{clara}}, \code{\link[stats]{kmeans}}.
#'    You can use \code{runSOM}, \code{runClara},
#'    \code{runPhenotype}, \code{runKmeans}, \code{runMclust} and
#'    \code{runHclust} to run clustering respectively.
#'
#' @export
#' @return An FSPY object with cluster
#'
#' @examples
#'
#' if (FALSE) {
#' # After building an FSPY object
#' # Set random seed to make results reproducible
#'
#' set.seed(1)
#' fspy <- runCluster(fspy, cluster.method = "som", xdim = 3, ydim = 3, verbose = TURE)
#'
#' # K-means clustering
#' fspy <- runCluster(fspy, cluster.method = "kmeans", k = 9, verbose = TRUE)
#'
#' # Clara clustering
#' fspy <- runCluster(fspy, cluster.method = "clara", k = 9, verbose = TRUE)
#'
#' # phenoGraph clustering
#' fspy <- runCluster(fspy, cluster.method = "phenograph", verbose = TRUE)
#'
#' # hclust clustering
#' # not recommended for large cell size
#' fspy <- runCluster(fspy, cluster.method = "hclust", k = 9, verbose = TRUE)
#'
#' # mclust clustering
#' # not recommended for large cell size
#' fspy <- runCluster(fspy, cluster.method = "mclust", verbose = TRUE)
#' }
#'
runCluster <- function(object, cluster.method = c("som", "kmeans", "clara", "phenograph", "hclust", "mclust"),
                       verbose = FALSE, ...) {

  if (missing(object)) {
    stop(Sys.time(), " [ERROR] FSPY object is missing ")
  }
  cluster.method <- match.arg(cluster.method)
  if (cluster.method == "som") {
    object <- runSOM(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$som.id
  } else if (cluster.method == "hclust") {
    object <- runHclust(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$hclust.id
  } else if (cluster.method == "mclust") {
    object <- runMclust(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$mclust.id
  } else if (cluster.method == "clara") {
    object <- runClara(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$clara.id
  } else if (cluster.method == "kmeans") {
    object <- runKmeans(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$kmeans.id
  } else if (cluster.method == "phenograph") {
    object <- runPhenograph(object, verbose = verbose, ...)
    object@meta.data$cluster.id <- object@meta.data$phenograph.id
  } else {
    warning(Sys.time(), " [WARNING] Invalid cluster.method parameter ")
  }

  # Initialization for root cells
  object@network <- list()
  object@meta.data$is.root.cells <- 0
  object@meta.data$is.leaf.cells <- 0

  return(object)

}

#'
#' processingCluster
#'
#' @name processingCluster
#'
#' @description Calculate Principal Components Analysis (PCA), t-Distributed
#'    Stochastic Neighbor Embedding (tSNE), Diffusion Map and Uniform Manifold
#'    Approximation and Projection (UMAP) of clusters calculated by runCluster.
#'
#' @param object an FSPY object
#' @param perplexity numeric. Perplexity parameter (should not be bigger than 3 *
#'    perplexity < nrow(X) - 1, see details for interpretation). See \code{\link[Rtsne]{Rtsne}}
#'    for more information.
#' @param k numeric. The parameter k in k-Nearest Neighbor.
#' @param downsampling.size numeric. Percentage of sample size of downsampling.
#'    This parameter is from 0 to 1. by default is 1.
#' @param force.resample logical. Whether to do resample if downsampling.size < 1
#' @param random.cluster logical. Whether to perfrom random downsampling. If FALSE, 
#'    an uniform downsampling will be processed.
#' @param umap.config object of class umap.config. See \code{\link[umap]{umap}}.
#' @param verbose logic. Whether to print calculation progress.
#' @param ... options to pass on to the dimensionality reduction functions.
#'
#' @seealso \code{\link[umap]{umap}}, \code{\link[gmodels]{fast.prcomp}},
#'    \code{\link[Rtsne]{Rtsne}}, \code{destiny}
#'
#' @return An FSPY object with cluster.id in meta.data
#'
#' @importFrom stats cutree
#'
#' @export
#' @return An FSPY object with dimensionality reduction of clusters
#'
#' @examples
#'
#' if (FALSE) {
#'
#' # After running clustering
#' set.seed(1)
#' fspy <- runCluster(fspy, cluster.method = "som", xdim = 3, ydim = 3, verbose = T)
#'
#' # Do not perfrom downsampling
#' fspy <- processingCluster(fspy, perplexity = 2)
#'
#' # Perform cluster based downsampling
#' # Only keep 50% cells
#' fspy <- processingCluster(fspy, perplexity = 2, downsampling.size = 0.5)
#'
#' # Processing clusters without downsampling step
#' fspy <- processingCluster(fspy, perplexity = 2, force.resample = FALSE)
#'
#' }
#'
processingCluster <- function(object, perplexity = 5, k = 5,
                              downsampling.size = 1,
                              force.resample = TRUE,
                              random.cluster = FALSE, 
                              umap.config = umap.defaults, verbose = FALSE,
                              ...) {

  if (missing(object)) {
    stop(Sys.time(), " [ERROR] FSPY object is missing ")
  }

  if (!"cluster.id" %in% colnames(object@meta.data)) {
    stop(Sys.time(), " [ERROR] cluster.id is not in colnames of FSPY object, please run runCluster first ")
  }

  # checking index of markers in cluster
  cluster.meta <- fetchClustMeta(object, verbose = FALSE)
  cluster.mat <- cluster.meta[, match(object@markers, colnames(cluster.meta))]

  # run PCA
  if (verbose) message(Sys.time(), " [INFO] Calculating PCA")
  pca.info <- fast.prcomp( t(cluster.mat), ...)
  colnames(pca.info$rotation) <- paste0("PC_", 1:ncol(pca.info$rotation))
  if (verbose) message(Sys.time(), " [INFO] Calculating tSNE")
  tsne.info <- Rtsne(as.matrix(cluster.mat), perplexity = perplexity, ...)
  colnames(tsne.info$Y) <- paste0("tSNE_", 1:ncol(tsne.info$Y))
  if (verbose) message(Sys.time(), " [INFO] Calculating Diffusion Map")
  dm.info <- DiffusionMap(cluster.mat, k=5, ...)
  colnames(dm.info@eigenvectors) <- paste0("DC_", 1:ncol(dm.info@eigenvectors))
  if (verbose) message(Sys.time(), " [INFO] Calculating UMAP")
  umap.config$n_neighbors <- k
  umap.info <- umap(cluster.mat, config = umap.config, ...)
  colnames(umap.info$layout) <- paste0("UMAP_", 1:ncol(umap.info$layout))

  object@cluster <- data.frame(pca.info$rotation, tsne.info$Y, dm.info@eigenvectors, umap.info$layout)
  rownames(object@cluster) <- rownames(object@tree.meta$cluster)

  if (force.resample) {
    # Initialization
    object@network <- list()
    object@meta.data$is.root.cells <- 0
    object@meta.data$is.leaf.cells <- 0
    object@meta.data$dowsample <- 0
    object@meta.data$pseudotime <- 0
    object@meta.data$traj.value <- 0
    object@meta.data$traj.value.log <- 0
    object@meta.data$is.root.cells <- 0
    object@meta.data$is.leaf.cells <- 0
    object@meta.data$branch.id <- 0
    object@pca.sdev <- vector()
    object@umap.value <- object@tsne.value <- object@pca.scores <- object@pca.value <- matrix()
    object@dm <- new("DiffusionMap")

    cell.sub <- NULL
    if (downsampling.size >= 1) {
      if (verbose) message(Sys.time(), " [INFO] No downsampling performed")
      cell.name <- object@meta.data$cell
    } else if ( downsampling.size <= 0) {
      warning(Sys.time(), " [WARNING] The value of downsampling.size must be larger than 0 ")
      cell.name <- object@meta.data$cell
    } else {
      if (random.cluster) {
        cell.name <- sapply(unique(object@meta.data$cluster.id), function(x) sample(object@meta.data$cell[which(object@meta.data$cluster.id == x)], ceiling(length(which(object@meta.data$cluster.id == x)) * downsampling.size )) )
        cell.name <- unlist(cell.name)
      } else {
        cell.name <- sapply(unique(object@meta.data$cluster.id), function(x) { 
          cell.sub <- as.character(object@meta.data$cell[which(object@meta.data$cluster.id == x)])
          cell.sub <- cell.sub[seq(1, length(cell.sub), by = 1/downsampling.size)]
          } )
        cell.name <- unlist(cell.name)
      }
      
    }

    object@cell.name <- as.character(cell.name)
    object@meta.data$dowsample[match(cell.name, object@meta.data$cell)] <- 1

  }

  return(object)

}


#'
#' runHclust
#'
#' @name runHclust
#'
#' @description Hierarchical cluster analysis on a set of dissimilarities
#'    and methods for analyzing it.
#'
#' @param object an FSPY object
#' @param hclust.method character or a function. The agglomeration method to be used.
#'    This should be one of "ward.D", "ward.D2", "single", "complete", "average",
#'    "mcquitty", "median" or "centroid". Or you can specify an equation as input, for example
#'    \code{function(x) hclust(x,method = 'ward.D2')}.
#' @param dist.method character or a function. The distance measure to be used.
#'    This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary"
#'    or "minkowski". Or you can specify an equation as input, for example
#'    \code{function(x) as.dist((1-cor(t(x)))/2)}.
#' @param k numeric. The number of clusters.
#' @param verbose logical. Whether to print calculation progress.
#'
#' @seealso \code{\link[stats]{hclust}}, \code{\link[stats]{dist}}
#'
#' @importFrom stats hclust dist
#'
#' @export
#' @return An FSPY object with cluster
#'
#' if (FALSE) {
#' fspy <- runHclust(fspy, k = 9, verbose = TRUE)
#' }
#'
#'
runHclust <- function(object, k = 25,
                      hclust.method = "complete", dist.method = "euclidean",
                      verbose = FALSE) {

  if (verbose) message(Sys.time(), " [INFO] Calculating Hclust.")

  # check dist parameters
  if (is.character(dist.method)) {
    d <- stats::dist(object@log.data, method = dist.method)
  } else if (is.function(dist.method)) {
    d <- dist.method(object@log.data)
  } else {
    warning(Sys.time(), " [WARNING] Invalid dist.method parameter.")
    d <- stats::dist(object@log.data)
  }

  # check hclust parameters
  if (is.character(hclust.method)) {
    hc <- stats::hclust(d, method = hclust.method)
  } else if (is.function()) {
    hc <- dist.method(d)
  } else {
    warning(Sys.time(), " [WARNING] Invalid hclust.method parameter.")
    hc <- stats::hclust(d)
  }


  hc.tree <- cutree(hc, k = k)

  object@meta.data$hclust.id <- object@meta.data$cluster.id <- hc.tree

  if (verbose) message(Sys.time(), " [INFO] Calculating Hclust completed.")
  return(object)
}



#'
#' runKmeans
#'
#' @name runKmeans
#'
#' @description Perform k-means clustering on a data matrix.
#'
#' @param object  an FSPY object
#' @param k numeric. The number of clusters.
#' @param iter.max numeric. The maximum number of iterations allowed.
#' @param nstart numeric. If k is a number, how many random sets should be chosen.
#' @param algorithm character. Type of algorithm that will be choosen to calculate 
#'    kmeans. Four algoritms are provided: Hartigan-Wong, Lloyd, Forgy, MacQueen.
#' @param trace logical or integer number.
#' @param scale logical. Whether to use scaled data in kmeans.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[stats]{kmeans}} function
#'
#' @return an FSPY object with kmeans.id in meta.data
#'
#' @seealso \code{\link[stats]{kmeans}}
#'
#' @importFrom stats kmeans
#' @export
#' @examples
#'
#' if (FALSE) {
#' fspy <- runKmeans(fspy, k = 25, verbose = TRUE)
#' }
#'
runKmeans <- function(object, k = 25, iter.max = 10, nstart = 1,
                      algorithm = c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"),
                      trace=FALSE, scale = FALSE, verbose = FALSE, ...) {

  if (verbose) message(Sys.time(), " [INFO] Calculating Kmeans.")

  if (scale) kmeans.data <- scale(object@log.data) else kmeans.data = object@log.data
  
  algorithm <- match.arg(algorithm)
  kmeans.info <- kmeans(kmeans.data, centers = k, iter.max = iter.max, nstart = nstart,
                        algorithm = algorithm, trace = FALSE)

  object@meta.data$kmeans.id <- object@meta.data$cluster.id  <- kmeans.info$cluster

  if (verbose) message(Sys.time(), " [INFO] Calculating Kmeans completed.")
  return(object)
}


#'
#' runClara
#'
#' @name runClara
#'
#' @description Clustering a data matrix into k clusters
#'
#' @param object  an FSPY object
#' @param k numeric. The number of clusters. It is required that
#'    0 < k < n where n is the number of observations (i.e., n = nrow(x)).
#' @param metric character. string specifying the metric to be used for
#'    calculating dissimilarities between observations.
#' @param stand logical. Indicating if the measurements in x are
#'    standardized before calculating the dissimilarities.
#' @param samples numeric. Say N, the number of samples to be drawn from the dataset.
#'    The default is N = 5,
#' @param trace numberic. Indicating a trace level for diagnostic output during the algorithm
#' @param scale logical. Whether to use scaled data in kmeans.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[cluster]{clara}} function
#'
#' @return an FSPY object with clara.id in meta.data
#'
#' @seealso \code{\link[cluster]{clara}}
#'
#' @importFrom cluster clara
#' @export
#' @examples
#' if (FALSE) {
#' fspy <- runClara(fspy, k = 25, verbose = TRUE)
#' }
#'
runClara <- function(object, k = 25, metric = c("euclidean", "manhattan", "jaccard"),
                     stand = FALSE, samples = 5, scale = TRUE,
                     trace = 0, verbose = FALSE, ...) {

  if (verbose) message(Sys.time(), " [INFO] Calculating Clara")

  if (scale) clara.data <- scale(object@log.data) else clara.data = object@log.data

  metric <- match.arg(metric)
  clara.info <- clara(clara.data, k = k, metric = metric, stand = stand, samples = samples,
                      trace = trace, ...)

  object@meta.data$clara.id <- object@meta.data$cluster.id  <- clara.info$clustering

  if (verbose) message(Sys.time(), " [INFO] Calculating Clara completed.")
  return(object)
}

#'
#' runMclust
#'
#' @name runMclust
#'
#' @description Model-based clustering based on parameterized finite Gaussian mixture models.
#'    This function is based on \code{\link[mclust]{Mclust}}.
#'
#' @param object  an FSPY object
#' @param scale logical. Whether to use scaled data in Mclust.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[mclust]{Mclust}} function
#'
#' @return an FSPY object with mclust.id in meta.data
#'
#' @seealso \code{\link[mclust]{Mclust}}
#'
#' @export
#'
#' @importFrom mclust Mclust mclustBIC
#' @examples
#' if (FALSE) {
#' fspy <- runMclust(fspy, verbose = TRUE)
#' }
#'
runMclust <- function(object, scale = FALSE,
                      verbose = FALSE, ...) {

  if (verbose) message(Sys.time(), " [INFO] Calculating Mclust.")

  if (scale) mclust.data <- scale(object@log.data) else mclust.data = object@log.data

  mod <- Mclust(mclust.data, ...)

  object@meta.data$mclust.id <- object@meta.data$cluster.id  <- mod$classification

  if (verbose) message(Sys.time(), " [INFO] Calculating Mclust completed.")
  return(object)
}



#'
#' calculation SOM in FSPY object
#'
#' @description Build a self-organizing map
#'
#' @param object  an FSPY object
#' @param xdim  Width of the grid.
#' @param ydim  Hight of the grid.
#' @param rlen  Number of times to loop over the training data for each MST
#' @param mst   Number of times to build an MST
#' @param alpha Start and end learning rate
#' @param radius Start and end radius
#' @param init  Initialize cluster centers in a non-random way
#' @param distf Distance function (1=manhattan, 2=euclidean, 3=chebyshev,
#'              4=cosine)
#' @param codes Cluster centers to start with
#' @param importance array with numeric values. Parameters will be scaled
#'                   according to importance
#' @param method the distance measure to be used. This must be one of "euclidean",
#'      "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'      Any unambiguous substring can be given. See \code{\link[stats]{dist}}
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[FlowSOM]{SOM}} function
#'
#' @return an FSPY object with som.id in FSPY object
#' @seealso \code{\link{BuildSOM}}
#'
#' @references This code is strongly based on the \code{\link[FlowSOM]{SOM}} function.
#'             Which is developed by Sofie Van Gassen, Britt Callebaut and Yvan Saeys (2018).
#'
#' @importFrom FlowSOM SOM
#'
#' @seealso \code{\link[FlowSOM]{SOM}}
#'
#' @export
#'
#' @examples
#' if (FALSE) {
#' fspy <- runSOM(fspy, xdim = 10, ydim = 10, verbose = TRUE)
#' }
#'
runSOM <- function(object, xdim = 6, ydim = 6, rlen = 8, mst = 1,
                   alpha = c(0.05,  0.01), radius = 1, init = FALSE,
                   distf = 2, codes = NULL, importance = NULL,
                   method = "euclidean", verbose= FALSE, ...) {

  if (verbose) message(Sys.time(), " [INFO] Calculating FlowSOM.")
  # FlowSOM
  flowset <- as.matrix(object@log.data)
  flowsom <- FlowSOM::SOM(flowset,
                          xdim = xdim, ydim = ydim, rlen = rlen, mst = mst,
                          alpha = alpha[1], radius = radius,
                          init = init,
                          distf = distf, silent = verbose,
                          codes = codes, importance = importance,
                          ...)

  # generating som network
  object@meta.data$som.id <- object@meta.data$cluster.id  <- flowsom$mapping[, 1]
  object@meta.data$som.value <- flowsom$mapping[, 2]
  object@som <- flowsom

  #object@som.network <- buildSOMnet(flowsom, object, method = method)

  if (verbose) message(Sys.time(), " [INFO] Calculating FlowSOM completed.")
  return(object)
}

#' RphenoGraph clustering
#'
#' @description
#'    A simple R implementation of the phenograph
#'    [PhenoGraph](http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm,
#'    which is a clustering method designed for high-dimensional single-cell
#'    data analysis. It works by creating a graph ("network") representing
#'    phenotypic similarities between cells by calculating the Jaccard
#'    coefficient between nearest-neighbor sets, and then identifying communities
#'    using the well known [Louvain method](https://sites.google.com/site/findcommunities/)
#'    in this graph.
#'
#' @param object an FSPY object.
#' @param scale logical. Whether to scale the expression matrix
#' @param knn numeric. Number of nearest neighbours, default is 30.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{igraph} function
#'
#'
#' @importFrom igraph graph.adjacency simplify distances
#' @return An FSPY object with cluster
#'
#' @export
#' @examples
#' if (FALSE) {
#' fspy <- runPhenograph(fspy, knn = 30, verbose = TRUE)
#' }
#'
runPhenograph <- function(object, knn = 30, scale = FALSE, verbose = FALSE, ...){


  if (verbose) message(Sys.time(), " [INFO] Calculating phenoGraph")

  if (scale) phenograph.data <- scale(object@log.data) else phenograph.data = object@log.data

  mod <- Rphenograph(phenograph.data, k = 30)

  object@meta.data$phenograph.id <- object@meta.data$cluster.id  <- membership(mod[[2]])

  if (verbose) message(Sys.time(), " [INFO] Calculating phenoGraph completed.")

  return(object)

}

#'
#' RphenoGraph clustering
#'
#' @name Rphenograph
#'
#' @description  R implementation of the PhenoGraph algorithm
#'
#' A simple R implementation of the [PhenoGraph](http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm,
#' which is a clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph ("network") representing
#' phenotypic similarities between cells by calclating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities
#' using the well known [Louvain method](https://sites.google.com/site/findcommunities/) in this graph.
#'
#' This function is developed by Hao Chen  and updated by Yuting Dai.
#'
#' @param data matrix; input data matrix
#' @param k integer; number of nearest neighbours (default:30)
#'
#' @return a list contains an igraph graph object for \code{graph_from_data_frame}
#'     and a communities object, the operations of this class contains:
#' \item{print}{returns the communities object itself, invisibly.}
#' \item{length}{returns an integer scalar.}
#' \item{sizes}{returns a numeric vector.}
#' \item{membership}{returns a numeric vector, one number for each vertex in
#'     the graph that was the input of the community detection.}
#' \item{modularity}{returns a numeric scalar.}
#' \item{algorithm}{returns a character scalar.}
#' \item{crossing}{returns a logical vector.}
#' \item{is_hierarchical}{returns a logical scalar.}
#' \item{merges}{returns a two-column numeric matrix.}
#' \item{cut_at}{returns a numeric vector, the membership vector of the vertices.}
#' \item{as.dendrogram}{returns a dendrogram object.}
#' \item{show_trace}{returns a character vector.}
#' \item{code_len}{returns a numeric scalar for communities found with the InfoMAP
#'     method and NULL for other methods.}
#' \item{plot}{for communities objects returns NULL, invisibly.}
#'
#' @references Jacob H. Levine and et.al. Data-Driven Phenotypic Dissection of AML
#'     Reveals Progenitor-like Cells that Correlate with Prognosis. Cell, 2015.
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' Rphenograph_out <- Rphenograph(data, k = 45)
#' modularity(Rphenograph_out[[2]])
#' membership(Rphenograph_out[[2]])
#' iris_unique$phenograph_cluster <- factor(membership(Rphenograph_out[[2]]))
#'
#' @importFrom igraph graph.data.frame cluster_louvain modularity membership
#' @import ggplot2
#' @import Rcpp
#' @useDynLib flowSpy
#'
#' @export
#' @return cluster information
#'
#' @author Hao Chen <chen_hao@immunol.a-star.edu.sg>
#'
Rphenograph <- function(data, k=30){
  if(is.data.frame(data))
    data <- as.matrix(data)

  if(!is.matrix(data))
    stop("Wrong input data, should be a data frame of matrix!")

  if(k<1){
    stop("k must be a positive integer!")
  }else if (k > nrow(data)-2){
    stop("k must be smaller than the total number of points!")
  }

  message("Run Rphenograph starts:","\n",
          "  -Input data of ", nrow(data)," rows and ", ncol(data), " columns","\n",
          "  -k is set to ", k)

  # cat("  Finding nearest neighbors...")
  t1 <- system.time(neighborMatrix <- find_neighbors(data, k=k+1)[,-1])
  # cat("DONE ~",t1[3],"s\n", " Compute jaccard coefficient between nearest-neighbor sets...")
  t2 <- system.time(links <- jaccard_coeff(neighborMatrix))

  # cat("DONE ~",t2[3],"s\n", " Build undirected graph from the weighted links...")
  links <- links[links[,1]>0, ]
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  t3 <- system.time(g <- graph.data.frame(relations, directed = FALSE))

  # Other community detection algorithms:
  #    cluster_walktrap, cluster_spinglass,
  #    cluster_leading_eigen, cluster_edge_betweenness,
  #    cluster_fast_greedy, cluster_label_prop
  # cat("DONE ~",t3[3],"s\n", " Run louvain clustering on the graph ...")
  t4 <- system.time(community <- cluster_louvain(g))
  # cat("DONE ~",t4[3],"s\n")

  message("Run Rphenograph DONE, totally takes ", sum(c(t1[3],t2[3],t3[3],t4[3])), "s.")
  # cat("  Return a community class\n  -Modularity value:", modularity(community),"\n")
  # cat("  -Number of clusters:", length(unique(membership(community))))

  return(list(g, community))
}


#' K Nearest Neighbour Search
#'
#' @name find_neighbors
#' @description  Uses a kd-tree to find the p number of
#'    near neighbours for each point in an input/output dataset.
#'   Use the nn2 function from the RANN package, utilizes the
#'   Approximate Near Neighbor (ANN) C++ library, which can give
#'   the exact near neighbours or (as the name suggests) approximate near neighbours
#'   to within a specified error bound. For more information on the ANN library please
#'   visit http://www.cs.umd.edu/~mount/ANN/.
#'
#' @param data matrix; input data matrix
#' @param k integer; number of nearest neighbours
#'
#' @return a n-by-k matrix of neighbor indices
#'
#' @author Hao Chen <chen_hao@immunol.a-star.edu.sg>
#'
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' neighbors <- find_neighbors(data, k=10)
#'
#' @importFrom RANN nn2
#' @export
#'
#'
find_neighbors <- function(data, k){
  nearest <- nn2(data, data, k, searchtype = "standard")
  return(nearest[[1]])
}







