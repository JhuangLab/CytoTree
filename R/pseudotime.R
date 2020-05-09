#'
#' definition of root cells
#'
#' @name defRootCells
#'
#' @description definition of root cells
#'
#' @param object an FSPY object
#' @param root.cells vector. Cell name of the root cells
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return An FSPY object
#'
#' @export
#'
#' @examples
#'
#' if (FALSE) {
#' # Define root cells by cluster
#' fspy <- defRootCells(fspy, root.cells = 6, verbose = TRUE)
#' fspy <- defRootCells(fspy, root.cells = c(6,8), verbose = TRUE)
#'
#' # Define root cells by cell names
#' cells <- test.meta.data$cell[which(test.meta.data$stage == "D0")]
#' cells <- as.character(cells)
#' fspy <- defRootCells(fspy, root.cells = cells, verbose = TRUE)
#' }
#'
#'
defRootCells <- function(object, root.cells = NULL, verbose = FALSE) {
  if (length(object@root.cells) != 0) message(Sys.time(), " [INFO] root.cells in FSPY object exist, they will be replaced.")

  if (!is.vector(root.cells)) stop(Sys.time(), " [ERROR] root.cells must be a vector")

  if (is.character(root.cells)) {
    root.cells <- root.cells[root.cells %in% object@meta.data$cell]
  } else if (is.numeric(root.cells)) {
    root.cells <- object@meta.data$cell[object@meta.data$cluster.id %in% root.cells]
  } else {
    stop(Sys.time(), " [ERROR] invalid root.cells .")
  }

  ds.cells <- object@meta.data$cell[which(object@meta.data$dowsample == 1)]
  root.cells <- root.cells[root.cells %in% ds.cells]

  object@meta.data$is.root.cells <- 0
  object@meta.data$is.root.cells[match(root.cells, object@meta.data$cell)] <- 1
  if ( length(root.cells) == 0 ) {
    stop(Sys.time(), " [ERROR] root.cells are not in meta.data")
  } else {
    object@root.cells <- root.cells
  }

  if (verbose) message(Sys.time(), " [INFO] ", length(root.cells),  " cells will be added to root.cells .")

  return(object)
}

#'
#' definition of leaf cells
#'
#' @name defLeafCells
#' @description definition of root cells
#'
#' @param object an FSPY object
#' @param leaf.cells character or numeric. Cell name of the root cells or
#'     cluster.id of root.cells
#' @param pseudotime.cutoff numeric. Cutoff of pseudotime. Cells with pseudotime
#'     over pseudotime.cutoff will be set to be leaf cells
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return An FSPY object
#'
#' @export
#'
#' @examples
#'
#' if (FALSE) {
#' # Define leaf cells by cluster
#' fspy <- defLeafCells(fspy, leaf.cells = 1, verbose = TRUE)
#' fspy <- defLeafCells(fspy, leaf.cells = c(1,3), verbose = TRUE)
#'
#' # Define root cells by cell names
#' cells <- test.meta.data$cell[which(test.meta.data$stage == "D10")]
#' cells <- as.character(cells)
#' fspy <- defLeafCells(fspy, leaf.cells = cells, verbose = TRUE)
#' }
#'
#'
defLeafCells <- function(object, leaf.cells = NULL, pseudotime.cutoff = 0, verbose = FALSE) {
  if (length(object@leaf.cells) != 0) message(Sys.time(), " [INFO] leaf.cells in FSPY object exist, they will be replaced.")

  if (!is.vector(leaf.cells)) stop(Sys.time(), " [ERROR] leaf.cells must be a vector")

  if (is.character(leaf.cells)) {
    leaf.cells <- leaf.cells[leaf.cells %in% object@meta.data$cell]
  } else if (is.numeric(leaf.cells)) {
    leaf.cells <- object@meta.data$cell[object@meta.data$cluster.id %in% leaf.cells]
  } else {
    stop(Sys.time(), " [ERROR] invalid leaf.cells.")
  }

  ds.cells <- object@meta.data$cell[which(object@meta.data$dowsample == 1)]
  leaf.cells <- leaf.cells[leaf.cells %in% ds.cells]

  if (pseudotime.cutoff > 0) {
    if ( !all("pseudotime" %in% colnames(object@meta.data)) ) {
      warning(Sys.time(), " [WARNING] pseudotime is not in meta.data of FSPY, please run Pseudotime first.")
      pseudotime.cutoff = 0
    }
  }

  leaf.time <- object@meta.data$pseudotime[match(leaf.cells, object@meta.data$cell)]
  leaf.cells <- leaf.cells[which(leaf.time >= pseudotime.cutoff )]

  object@meta.data$is.leaf.cells <- 0
  object@meta.data$is.leaf.cells[match(leaf.cells, object@meta.data$cell)] <- 1
  if ( length(leaf.cells) == 0 ) {
    stop(Sys.time(), " [ERROR] leaf.cells are not in meta.data")
  } else {
    object@leaf.cells <- leaf.cells
  }

  if (verbose) message(Sys.time(), " [INFO] ", length(leaf.cells),  " cells will be added to leaf.cells .")

  return(object)
}



#'
#' Calculation of Pseudotime
#'
#' @name runPseudotime
#'
#' @description calculation of Pseudotime based on KNN
#'
#' @param object An FSPY object
#' @param mode character. Specifies how igraph should interpret the supplied matrix.
#'    Possible values are: directed, undirected, upper, lower, max, min, plus.
#' @param dim.type character. Type of dimensionality reduction method used to calculate
#'    pseudotime: raw, umap, tsne, dc and pca. By default is raw.
#' @param dim.use numeric. Dimensions used to calculate pseudotime
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to calculation function.
#'
#' @importFrom igraph graph.adjacency simplify distances
#' @return An FSPY object
#'
#' @export
#' 
#' @importFrom stats kmeans
#'
#' @examples
#'
#' if (FALSE) {
#' 
#' fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "raw")
#' fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "umap", dim.use = 1:2)
#' fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "tsne", dim.use = 1:2)
#' fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "dc", dim.use = 1:3)
#' fspy <- runPseudotime(fspy, verbose = TRUE, dim.type = "pca", dim.use = 1:3)
#'
#' # tSNE plot colored by pseudotime
#' plot2D(fspy, item.use = c("tSNE_1", "tSNE_2"), category = "numeric",
#'        size = 1, color.by = "pseudotime") +
#'        scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))
#' # UMAP plot colored by pseudotime
#' plot2D(fspy, item.use = c("UMAP_1", "UMAP_2"), category = "numeric",
#'        size = 1, color.by = "pseudotime") +
#'        scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))
#' }
#'
runPseudotime <- function(object, mode = "undirected",
                          dim.type = c("raw", "pca", "tsne", "dc", "umap"), dim.use = 1:2,
                          verbose = FALSE, ...) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing.")

  if (verbose) message(Sys.time(), " [INFO] Calculating Pseudotime.")

  if ("pseudotime" %in% colnames(object@meta.data)) message(Sys.time(), " [INFO] Pseudotime exists in meta.data, it will be replaced.")

  if (missing(object)) stop(Sys.time(), " [ERROR] FSPY object is missing.")

  dim.type <- match.arg(dim.type)
  if (dim.type %in% c("tsne", "tSNE", "TSNE", "t-SNE","t_SNE", "t") ) {
    dim.name <- paste0("tSNE_", dim.use)
    mat <- object@tsne.value[, dim.name]
  } else if ( dim.type %in% c("PCA", "pca", "p") ) {
    dim.name <- paste0("PC_", dim.use)
    mat <- object@pca.value[, dim.name]
  } else if (dim.type %in% c("dc", "diffusionmap", "diffusion-map", "destiny", "d")) {
    dim.name <- paste0("DC_", dim.use)
    mat <- object@dm@eigenvectors[, dim.name]
  } else if (dim.type %in% c("umap", "UMAP", "u")) {
    dim.name <- paste0("UMAP_", dim.use)
    mat <- object@umap.value[, dim.name]
  } else {
    if (verbose) message(Sys.time(), " [INFO] The log data will be used to calculate pseudotime")
    mat <- object@log.data[which(object@meta.data$dowsample == 1), ]
  }

  object@meta.data$seed.pseudotime <- 0
  object@meta.data$core.pseudotime <- 0
  if (dim(mat)[1] > 40000) {
    kmean.mat <- kmeans(mat, centers = 40000)
    sub.cell <- kmean.mat$cluster[match(1:40000, kmean.mat$cluster)]
    object@meta.data$seed.pseudotime[match(names(sub.cell), object@meta.data$cell)] <- 1
    object@meta.data$core.pseudotime[match(names(kmean.mat$cluster), object@meta.data$cell)] <- kmean.mat$cluster
    sub.cell.name <- object@meta.data$cell[which(object@meta.data$seed.pseudotime == 1)]
    mat <- mat[match(sub.cell.name, rownames(mat)), ]
    object <- runKNN(object, given.mat = mat, verbose = FALSE)
  } else {
    object@meta.data$seed.pseudotime[match(rownames(mat), object@meta.data$cell)] <- 1
    object@meta.data$core.pseudotime[match(rownames(mat), object@meta.data$cell)] <- 1:dim(mat)[1]
    object <- runKNN(object, given.mat = mat, verbose = FALSE)
  }

  knn.index <- object@knn.index
  adj <- matrix(0, nrow(knn.index), nrow(knn.index))
  rownames(adj) <- colnames(adj) <- rownames(knn.index)
  for(i in seq_len(nrow(knn.index))) {
    adj[i, rownames(knn.index)[knn.index[i,]]] <- 1
  }


  g <- igraph::graph.adjacency(adj, mode = mode, ... )
  # remove self loops
  g <- simplify(g)

  root.cells <- object@root.cells[object@root.cells %in% rownames(knn.index)]
  dist.all.path <- distances(g, v = root.cells)
  dist.all.path[which(is.infinite(dist.all.path))] <- NA
  pst <- colMeans(dist.all.path, na.rm = TRUE)
  idx <- which(!is.na(pst))
  pst[idx] <- ( pst[idx] - min(pst[idx]) )/ max( pst[idx] - min(pst[idx]) )
  pst.1 <- pst
  names(pst.1) <- object@meta.data$core.pseudotime[match(names(pst), object@meta.data$cell)]
  pst.info <- object@meta.data$core.pseudotime[which(object@meta.data$core.pseudotime > 0)]
  pseudotime <- pst.1[match(pst.info, names(pst.1))]
    
  object@meta.data$pseudotime <- 0
  object@meta.data$pseudotime[which(object@meta.data$dowsample == 1)] <- pseudotime
  object@meta.data$traj.value <- 0
  object@meta.data$traj.value.log <- 0

  if (verbose) message(Sys.time(), " [INFO] Calculating Pseudotime completed.")

  return(object)
}









