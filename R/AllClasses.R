#' @import Matrix
#' @import ggplot2
#' @import gmodels
#' @import Rtsne
#' @import destiny
#' @import FlowSOM
#' @import BiocNeighbors
#' @import matrixStats
#' @import flowUtils
#' @import umap
#' @import prettydoc
#' @import scatterpie
NULL


#'
#' Class \code{CYT}
#'
#' @aliases CYTclass, CYT-class, CYT
#'
#' @description  All information stored in CYT object.
#'    You can use \code{creatCYT} to   create an CYT
#'    object. In this package, most of the functions will use
#'    CYT object as input, and return a modified CYT obejct as well.
#'
#' @slot raw.data matrix. Raw signal data captured in flow
#'     or mass cytometry.
#' @slot log.data matrix. Log-transfromed dataset of raw.data.
#' @slot meta.data data.frame. Meta data information, and
#'     colnames of "stage" and "cell" are required.
#' @slot markers vector. Markers used in the calculation
#'     of PCA, tSNE, diffusion map and UMAP.
#' @slot markers.idx vector. Index of markers used in the
#'     calculation of PCA, tSNE, destiny and umap.
#' @slot cell.name vector. Cell names after performing downsampling.
#' @slot knn numeric. Numbers of nearest neighbors
#' @slot knn.index,knn.distance matrix. Each row of the
#'     \code{knn.index} matrix corresponds to a point
#'     in \code{log.data} and contains the row indices in
#'     \code{log.data} that are its nearest neighbors.
#'     And each row of the \code{knn.distance} contains
#'     the distance of its nearest neighbors.
#' @slot som list. Store som network information calculated
#'     using \code{\link[FlowSOM]{FlowSOM}}.
#' @slot cluster data.frame. Cluster information
#' @slot pca.sdev,pca.value,pca.scores PCA information of CYT
#'     object which are generated from \code{\link[gmodels]{fast.prcomp}}.
#' @slot tsne.value matrix. tSNE coordinates information.
#'     See \code{\link[Rtsne]{Rtsne}}.
#' @slot dm DiffusionMap object. Diffusion map calculated by package 
#'     \code{destiny}
#' @slot umap.value matrix umap coordinates information
#'     calculated using \code{\link[umap]{umap}}.
#' @slot root.cells vector, Names of root cells, which can
#'     be modified by \code{defRootCells}.
#'     An root cell is manually set to be the origin of all cells.
#'     Pseudotime in root cells are the lowest.
#' @slot leaf.cells vector. Names of leaf cells, which can be
#'     modified by \code{defLeafCells}.
#'     An leaf cell is manually set to be the terminal state of
#'     all cells. Pseuodtime in leaf cells are the largest.
#' @slot network list. Network stored in the calculation of
#'     trajectory and pseudotime.
#' @slot walk list. Random forward and backward walk between
#'     \code{root.cells} and \code{leaf.cells}.
#' @slot diff.traj list. Differentiation trajectory all cells.
#' @slot plot.meta data.frame. Plot meta information for
#'     \code{plot2D} or \code{plot3D}.
#' @slot tree.meta data.frame. Tree meta information of CYT object.
#'
#' @importClassesFrom destiny DiffusionMap DPT
#'
#' @useDynLib CytoTree
#'
#' @export
#'
#'
#' @return NULL
#'
#'
setClass("CYT", slots = c(
    raw.data = "matrix",
    log.data = "matrix",
    meta.data = "data.frame",
    markers = "vector",
    markers.idx = "vector",
    cell.name = "vector",

    # KNN
    knn = "numeric",
    knn.index = "matrix",
    knn.distance = "matrix",

    # som network
    som = "list",

    # cluster information
    cluster = "data.frame",

    # pca information
    pca.sdev = "vector",
    pca.value = "matrix",
    pca.scores = "matrix",

    # tsne information
    tsne.value = "matrix",

    # diffusion map information
    dm = c("DiffusionMap", NULL),

    # umap information
    umap.value = "matrix",

    # run for pseudotime
    root.cells = "vector",
    leaf.cells = "vector",
    network = "list",

    # trajectory analysis
    walk = "list",
    diff.traj = "list",

    # for visualization
    plot.meta = "data.frame",
    tree.meta = "data.frame"
    )
)

setValidity("CYT", function(object) {
  X <- object@log.data
  n <- nrow(X)
  p <- ncol(X)
  if(!is.matrix(X)) {
    return("Log data must be matrix.")
  }
  return(TRUE)
})

#' create an CYT object
#'
#' @description This function is about how to build an CYT object.
#'    An CYT object is the base for the whole analysizing workflow
#'    of flow and mass cytometry data.
#'
#' @name createCYT
#'
#' @param raw.data matrix. Raw data read from FCS file after perform
#'    preprocessing.
#' @param markers vector. Detailed marker information in the gate of
#'    flow cytometer.
#' @param meta.data data.frame. Raw metadata of each cell.
#'    Columns "cell" and "stage" are required.
#' @param batch vector. Batch covariate (only one batch allowed).
#'    Method to correct batch effect
#'    function is refered to \code{\link[sva]{ComBat}}.
#' @param batch.correct logical. Whether to correct batch effect.
#'    If TRUE, batch must be provided.
#' @param normalization.method character. Normalization and transformation
#'    method. Whether to normalize and log transformed of raw.data.
#'    In CytoTree workflow, it's better to perform transformation of
#'    FCS data using \code{runExprsExtract} or \code{runExprsMerge}
#'    before creating an CYT object. \code{CytoTree} only provide
#'    log transforma method. If you need to using truncateTransform,
#'    scaleTransform, linearTransform, quadraticTransform and
#'    lnTransform, see \code{flowCore} for more
#'    information. And \code{runExprsExtract} in
#'    \code{CytoTree}, autoLgcl, cytofAsinh, logicle, arcsinh,
#'    and logAbs can be used to perform transformation of FCS data.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... paramters pass to \code{correctBatchCYT} function.
#'
#' @importFrom methods new
#' @importFrom stats median
#' @useDynLib CytoTree
#'
#' @export
#'
#' @return An CYT object with raw.data and markers and meta.data
#'
#' @examples
#'
#'
#' # Read fcs files
#' fcs.path <- system.file("extdata", package = "CytoTree")
#' fcs.files <- list.files(fcs.path, pattern = '.FCS$', full = TRUE)
#'
#' fcs.data <- runExprsMerge(fcs.files, comp = FALSE, transformMethod = "none")
#'
#' # Refine colnames of fcs data
#' recol <- c(`FITC-A<CD43>` = "CD43", `APC-A<CD34>` = "CD34", 
#'            `BV421-A<CD90>` = "CD90", `BV510-A<CD45RA>` = "CD45RA", 
#'            `BV605-A<CD31>` = "CD31", `BV650-A<CD49f>` = "CD49f",
#'            `BV 735-A<CD73>` = "CD73", `BV786-A<CD45>` = "CD45", 
#'            `PE-A<FLK1>` = "FLK1", `PE-Cy7-A<CD38>` = "CD38")
#' colnames(fcs.data)[match(names(recol), colnames(fcs.data))] = recol
#' fcs.data <- fcs.data[, recol]
#' 
#' day.list <- c("D0", "D2", "D4", "D6", "D8", "D10")
#' meta.data <- data.frame(cell = rownames(fcs.data),
#'                         stage = gsub(".FCS.+", "", rownames(fcs.data) ) )
#' meta.data$stage <- factor(as.character(meta.data$stage), levels = day.list)
#' 
#' markers <- c("CD43","CD34","CD90","CD45RA","CD31","CD49f","CD73","CD45","FLK1","CD38")
#' 
#'# Build the CYT object
#' cyt <- createCYT(raw.data = fcs.data, markers = markers,
#'                  meta.data = meta.data,
#'                  normalization.method = "log",
#'                  verbose = TRUE)
#' 
#' # See information
#' cyt
#'
#'
createCYT <- function(raw.data, markers, meta.data,
                      batch = NULL, batch.correct = FALSE,
                      normalization.method = "none",
                      verbose = FALSE, ...) {
  # QC of cells
  if (missing(raw.data)) stop(Sys.time(), " [ERROR] raw.data is required")
  if (!is.matrix(raw.data)) {
    warning(Sys.time(), " [WARNING] raw.data must be a matrix")
    raw.data <- as.matrix(raw.data)
  }
  if (verbose) message(Sys.time(), " [INFO] Number of cells in processing: ", dim(raw.data)[1])

  # QC of metadata
  if (missing(meta.data)) stop(Sys.time(), " [ERROR] meta.data must be a data.frame")
  if (!is.data.frame(meta.data)) {
    warning(Sys.time(), " [WARNING] meta.data must be a data.frame")
    meta.data <- as.matrix(meta.data)
  }

  if (!all(c("cell", "stage") %in% colnames(meta.data))) {
    stop(Sys.time(), " [ERROR] cell and stage information must be provided in meta.data")
  }

  if (nrow(raw.data) != nrow(meta.data)) {
    stop(Sys.time(), " [ERROR] cell number in raw.data is not equal to that in meta.data")
  } else {
    if (verbose) message(Sys.time(), " [INFO] rownames of meta.data and raw.data will be named using column cell")
    rownames(raw.data) = as.character(meta.data$cell)
    rownames(meta.data) = as.character(meta.data$cell)
  }

  # load index of markers of FCS
  if (missing(markers)) stop(Sys.time(), " [ERROR] markers is missing")
  if (!is.vector(markers)) {
    warning(Sys.time(), " [WARNING] markers must be a vector")
    markers <- as.vector(markers)
  }

  # check markers' index in raw.data
  markers.idx <- match(markers, colnames(raw.data))
  if (verbose) message(Sys.time(), " [INFO] Index of markers in processing")
  if (any(is.na(markers.idx))) {
    sub.markers <- markers[which(is.na(markers.idx))]
    warning(Sys.time(), " [WARNING] ", sub.markers, " not existes in colnames
            of raw.data. It will be removed. ")

    markers <- markers[which(!is.na(markers.idx))]
    markers.idx <- markers.idx[which(!is.na(markers.idx))]
  }

  # Create an CYT object
  if (verbose) message(Sys.time(), " [INFO] Creating CYT object.")
  object <- methods::new("CYT", raw.data = raw.data, meta.data = meta.data,
                markers = markers, markers.idx = markers.idx)

  # normalization and Log-normalize the data
  if (normalization.method == "log") {
    if (verbose) message(paste0(Sys.time(), " [INFO] Determining normalization factors"))
    all.log.data <- abs(raw.data)
    cs <- apply(all.log.data, 2, sum)
    norm_factors <- (10**ceiling(log10(median(cs))))/cs
    norm_factors_idx <- which(!is.na(norm_factors))
    if (length(which(is.na(norm_factors))) > 0) {
      warning(paste0(Sys.time(), " [WARNING] Unavailable log data column, please check your data"))
    }
    if (verbose) message(paste0(Sys.time(), " [INFO] Normalization and log-transformation."))
    object@raw.data[, norm_factors_idx] <- round(
      log10(sweep(all.log.data[, norm_factors_idx], 2,
                  norm_factors[norm_factors_idx], "*")+1), digits=3)

    object@log.data <- object@raw.data[, markers.idx]
  } else if (normalization.method == "none") {
    if (verbose) message(Sys.time(), " [INFO] No normalization and transformation ")
    object@log.data <- raw.data[, markers.idx]
  } else {
    if (verbose) message(Sys.time(), " [INFO] No normalization and transformation ")
    object@log.data <- raw.data[, markers.idx]
  }

  # correcting batch effect
  if (batch.correct) {
    if (is.null(batch)) {
      warning(Sys.time(), " [WARNING] batch must be provided when batch.correct is TRUE ")
    } else {
      object <- correctBatchCYT(object, batch = batch, ...)
    }
  }

  # Initialization of all parameters in computation
  object@plot.meta <- data.frame(row.names = object@meta.data$cell)
  object@meta.data$dowsample <- 1
  object@meta.data$pseudotime <- 0
  object@meta.data$traj.value <- 0
  object@meta.data$traj.value.log <- 0
  object@meta.data$is.root.cells <- 0
  object@meta.data$is.leaf.cells <- 0

  if (verbose) message(Sys.time(), " [INFO] Build CYT object succeed ")
  return(object)
}













