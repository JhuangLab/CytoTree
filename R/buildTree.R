#'
#' buildTree
#'
#' @name buildTree
#'
#' @param object a CYT object
#' @param method character. Mehtod to build MST.
#' @param dim.type character. Type of dimensions that will be used to build the tree.
#'    Five \code{dim.type} are provided, 'raw', 'pca', 'tsne', 'dc' and 'umap'. 
#'    By default is 'raw'.
#' @param dim.use numeric. Number of dimensions that will be used to build the tree.
#'    For example. If \code{dim.use} is 'raw', there is no limit for \code{dim.type}. 
#'    And if the \code{dim.use} is 'tsne' or 'umap', the default \code{dim.use} is seq_len(2).
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#'
#' @importFrom stats aggregate
#' @importFrom igraph cluster_louvain membership graph.adjacency minimum.spanning.tree
#'
#' @return A CYT object with tree
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' 
#' cyt <- buildTree(cyt, dim.type = "raw")
#'
#' # build minimum spanning tree (MST) based on tsne
#' cyt <- buildTree(cyt, dim.type = "tsne", dim.use = seq_len(2))
#'
#' # Using PCA
#' cyt <- buildTree(cyt, dim.type = "pca", dim.use =seq_len(4))
#'
#' # Using UMAP
#' cyt <- buildTree(cyt, dim.type = "umap", dim.use = seq_len(2))
#'
#' # Using Diffusion Maps
#' cyt <- buildTree(cyt, dim.type = "dc", dim.use = seq_len(3))
#' 
#'
buildTree <- function(object, method = "euclidean",
                      dim.type = c("raw", "pca", "tsne", "dc", "umap"), 
                      dim.use = seq_len(2),
                      verbose = FALSE) {

  if (verbose) message(Sys.time(), " Calculating buildTree.")
  if (missing(object)) stop(Sys.time(), " CYT object is missing.")

  dim.type <- match.arg(dim.type)
  if (dim.type %in% c("tsne", "tSNE", "TSNE", "t-SNE","t_SNE", "t") ) {
    dim.name <- paste0("tSNE_", dim.use)
    tree.mat <- object@tsne.value[, dim.name]
  } else if ( dim.type %in% c("PCA", "pca", "p") ) {
    dim.name <- paste0("PC_", dim.use)
    tree.mat <- object@pca.value[, dim.name]
  } else if (dim.type %in% c("dc", "diffusionmap", "diffusion-map", "destiny", "d")) {
    dim.name <- paste0("DC_", dim.use)
    tree.mat <- object@dm@eigenvectors[, dim.name]
  } else if (dim.type %in% c("umap", "UMAP", "u")) {
    dim.name <- paste0("UMAP_", dim.use)
    tree.mat <- object@umap.value[, dim.name]
  } else {
    if (verbose) message(Sys.time(), " The log data will be used to calculate trajectory")
    tree.mat <- object@log.data[which(object@meta.data$dowsample == 1), object@markers.idx]
  }

  if (! "cluster.id" %in% colnames(object@meta.data)) {
    stop(Sys.time(), " Invalid cluster.id, please run runCluster first")
  }
  cluster.info <- object@meta.data$cluster.id[which(object@meta.data$dowsample == 1)]
  mst.mat <- aggregate(tree.mat, list(cluster = cluster.info), mean)
  rownames(mst.mat) <- mst.mat[, 1]

  adjacency <- stats::dist(mst.mat[, -1], method = method)
  fullGraph <- igraph::graph.adjacency(as.matrix(adjacency),
                                       mode = "undirected",
                                       weighted = TRUE)
  fullGraph <- simplify(fullGraph)
  tree.graph <- igraph::minimum.spanning.tree(fullGraph)


  # identify branch
  object@meta.data$branch.id <- membership(cluster_louvain(tree.graph))[match(object@meta.data$cluster.id,names(membership(cluster_louvain(tree.graph))))]

  # Storing network information
  object@network <- list(mst = tree.graph,
                         method = method,
                         dim.type = dim.type,
                         dim.use = dim.use,
                         mst.mat = mst.mat,
                         branch.id = membership(cluster_louvain(tree.graph)))

  # Initialization for root.cells and leaf cells
  if (verbose) message(Sys.time(), " Initialization for root.cells and leaf cells")
  object@meta.data$is.root.cells <- 0
  object@meta.data$is.leaf.cells <- 0

  # update tree meta information
  object <- updateClustMeta(object, verbose = FALSE)

  if (verbose) message(Sys.time(), " Calculating buildTree completed.")
  return(object)
}











