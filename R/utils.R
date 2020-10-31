#'
#' Update plot meta information of CYT
#'
#' @name updatePlotMeta
#'
#' @param object A CYT object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#' @return A CYT object
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' 
#' cyt <- updatePlotMeta(cyt)
#' 
#'
#'
updatePlotMeta <- function(object, verbose = TRUE) {
  plot.meta <- object@meta.data[which(object@meta.data$dowsample == 1), ]
  if (dim(object@pca.value)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@pca.value)
  }
  if (dim(object@tsne.value)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@tsne.value)
  }
  if (dim(object@dm@eigenvectors)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@dm@eigenvectors)
  }
  if (dim(object@umap.value)[2] > 0) {
    plot.meta <- cbind(plot.meta, object@umap.value)
  }

  plot.meta <- as.data.frame(plot.meta)
  object@plot.meta <- plot.meta

  if (verbose) message(Sys.time(), " Columns can be used in plot2D and plot3D: ", paste(colnames(plot.meta), collapse = " "))

  return(object)
}

#' 
#' Change marker used in the calculation of CYT
#' 
#' @name changeMarker
#' 
#' @param object A CYT object
#' @param markers vector. Markers used in the calculation
#'     of PCA, tSNE, diffusion map and UMAP.
#' @param verbose logical. Whether to print calculation progress.
#'     
#' @export
#' @return A CYT object
#' 
#' @examples
#' 
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' 
#' markers <- c("CD43", "CD34", "CD90", "CD45RA")
#' 
#' cyt <- changeMarker(cyt, markers = markers)
#' 
#' 
#'       
changeMarker <- function(object, markers = NULL, verbose = FALSE) {
  
  if (is.null(markers)) {
    stop(Sys.time(), " Markers not input")
  }
  if (!is.vector(markers)) {
    warning(Sys.time(), " Markers must be a vector")
    markers <- as.vector(markers)
  }
  if ( (sum(!markers %in% object@markers) == 0) & (sum(!object@markers %in% markers) == 0) ) {
    if (verbose) message(Sys.time(), " Markers is not changed")
    return(object)
  } 
  
  # check markers' index in raw.data
  markers.idx <- match(markers, colnames(object@raw.data))
  if (verbose) message(Sys.time(), " Index of markers in processing")
  if (any(is.na(markers.idx))) {
    sub.markers <- markers[which(is.na(markers.idx))]
    warning(Sys.time(), " ", sub.markers, " not exist in colnames
            of raw.data. It will be removed. ")
    
    markers <- markers[which(!is.na(markers.idx))]
    markers.idx <- markers.idx[which(!is.na(markers.idx))]
  }
  object@markers <- markers
  object@markers.idx <- markers.idx
  
  return(object)
  
}

#' 
#' Add meta information of CYT
#' 
#' @name addMetaData
#' 
#' @param object A CYT object
#' @param meta.info a vector, meta data of cell information
#' @param name character, colname of `meta.info`
#' @param verbose logical. Whether to print calculation progress.
#' 
#' @export
#' @return A CYT object
#' 
#' @examples
#' 
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' 
#' plot.meta <- fetchPlotMeta(cyt)
#' meta.info <- 1:nrow(plot.meta)
#' names(meta.info) <- plot.meta$cell
#' 
#' cyt <- addMetaData(cyt, meta.info = meta.info, name = "MyInformation")
#' plot.meta <- fetchPlotMeta(cyt)
#' 
#' 
addMetaData <- function(object, meta.info,
                        name = "NewCol",
                        verbose = FALSE) {
  
  if (missing(meta.info)) {
    stop(Sys.time(), " meta.info is missing")
  } 
  if (length(meta.info) != nrow(object@raw.data)) {
    stop(Sys.time(), " meta.info must be a vector or factor and the length is equal to cell number")
  }
  if (is.null(names(meta.info))) {
    if (verbose) message(Sys.time(), " the name of meta.info is missing")
    names(meta.info) <- rownames(object@meta.data)
  }
  if (name %in% colnames(object@meta.data)) {
    if (verbose) message(Sys.time(), " the colname is exist in meta.data of CYT object. The old one will be replaced")
    object@meta.data[, which(colnames(object@meta.data) == name)] <- meta.info[match(rownames(object@meta.data), names(meta.info))]
  } else {
    object@meta.data$NewCol <- meta.info[match(rownames(object@meta.data), names(meta.info))]
    colnames(object@meta.data)[which(colnames(object@meta.data) == "NewCol")] = name
  }
  
  return(object)
}




#'
#' Update clusters' meta information of CYT
#'
#' @name updateClustMeta
#'
#' @param object A CYT object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#' @return A CYT object
#'
#' @importFrom stats aggregate
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' 
#' cyt <- updateClustMeta(cyt)
#' 
#'
updateClustMeta <- function(object, verbose = FALSE) {

  # Generating tree meta information
  plot.data <- fetchPlotMeta(object, verbose = FALSE)
  plot.data <- cbind(plot.data, object@log.data[which(object@meta.data$dowsample == 1), ])

  if (length(unique(plot.data$stage)) > 1) {
    cell.count <- table(plot.data[, match(c("cluster.id", "stage"), colnames(plot.data)) ])
    cell.count<- t(sapply(seq_len(dim(cell.count)[1]), function(x) cell.count[x, ]))
    cell.total.number <- rowSums(cell.count)
    cell.total.number.percent <- cell.total.number/sum(cell.total.number)
    cell.percent<- t(sapply(seq_len(dim(cell.count)[1]), function(x) cell.count[x,]/sum(cell.count[x,]) ))
    colnames(cell.percent) <- paste0(colnames(cell.count), ".percent")
  } else {
    cell.count <- table(plot.data[, match(c("cluster.id", "stage"), colnames(plot.data)) ])
    cell.count<- t(sapply(seq_len(dim(cell.count)[1]), function(x) cell.count[x, ]))
    cell.count <- as.data.frame(t(cell.count))
    cell.total.number <- rowSums(cell.count)
    cell.total.number.percent <- cell.total.number/sum(cell.total.number)
    cell.percent<- t(sapply(seq_len(dim(cell.count)[1]), function(x) cell.count[x,]/sum(cell.count[x,]) ))
    cell.percent <- as.data.frame(t(cell.percent))
    colnames(cell.count) <- paste0(unique(plot.data$stage))
    colnames(cell.percent) <- paste0(unique(plot.data$stage), ".percent")
  }

  idx.redim <- match(c(colnames(object@log.data), "pseudotime", "traj.value", "traj.value.log"), colnames(plot.data))
  idx.redim <- unique(idx.redim)
  tree.meta <- stats::aggregate(plot.data[, idx.redim], list(cluster = plot.data[, "cluster.id"]), mean)
  tree.meta.1 <- data.frame(cell.count,
                           cell.number = cell.total.number,
                           cell.number.percent = cell.total.number.percent,
                           cell.percent)
  tree.meta <- cbind(tree.meta, tree.meta.1)
  tree.meta$branch.id <- plot.data$branch.id[match(tree.meta$cluster, plot.data$cluster.id)]

  object@tree.meta <- tree.meta
  if (verbose) message(Sys.time(), " Columns can be used in plotTree: ", paste(colnames(tree.meta), collapse = " "))

  return(object)
}


#'
#' Fetching plot metadata of CYT
#'
#' @name fetchPlotMeta
#'
#' @param object An CYT object
#' @param markers vector. Makers fetched from expression matrix
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return a data.frame containing meta information for visualization
#'
#' @export
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' 
#' plot.data <- fetchPlotMeta(cyt)
#' head(plot.data)
#'
#' plot.data <- fetchPlotMeta(cyt, markers = c("CD43", "CD34"))
#' head(plot.data)
#' 
#'
#'
#'
fetchPlotMeta <- function(object, markers = NULL, verbose = FALSE) {

    # update and fetch plot meta information
  object <- updatePlotMeta(object, verbose = FALSE)
  plot.meta <- object@plot.meta
  idx <- match(markers, colnames(object@log.data))
  idx <- idx[which(!is.na(idx))]
  if (length(idx) > 0) {
    plot.meta <- cbind(plot.meta, object@log.data[which(object@meta.data$dowsample == 1), idx])
  }


  return(plot.meta)
}

#'
#' Fetching clusters' metadata of CYT
#'
#' @name fetchClustMeta
#'
#' @param object An CYT object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return a data.frame containing clustering information for visualization
#'
#' @export
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' 
#' clust.data <- fetchClustMeta(cyt)
#' head(clust.data)
#' 
#'
#'
fetchClustMeta <- function(object, verbose = FALSE) {

  object <- updateClustMeta(object, verbose = verbose)

  return(object@tree.meta)
}


#'
#' Fetching cellls of CYT
#'
#' @name fetchCell
#'
#' @param object An CYT object
#' @param logical.connect character. "and" or "or"
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Paramters to pass to limitation
#'
#' @return a vector containing cell names
#'
#' @export
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' 
#' cell.fetch <- fetchCell(cyt, traj.value.log = 0.01)
#' cell.fetch <- fetchCell(cyt, stage = c("D0", "D10"))
#' cell.fetch <- fetchCell(cyt, stage = c("D0", "D10"), traj.value.log = 0.01,
#'                         logical.connect = "or")
#' 
#'
fetchCell <- function(object, logical.connect = "or", verbose = FALSE, ... ) {

  object <- updatePlotMeta(object, verbose = verbose)

  param.list <- list(...)
  plot.meta <- object@plot.meta

  cell.left <- NULL

  if (length(param.list) > 0) {
    param.list <- param.list[names(param.list) %in% colnames(plot.meta)]
    for (i in seq_along(param.list)) {
      sub <- param.list[[i]]
      sub.name <- names(param.list)[i]
      if (sub.name %in% c("cell", "stage")) {
        cell.sub <- plot.meta[plot.meta[, sub.name] %in% sub, "cell"]
      } else if (grepl(".id$", sub.name)) {
        cell.sub <- plot.meta[plot.meta[, sub.name] %in% sub, "cell"]
      } else if (is.numeric(sub)) {
        if (length(sub) == 1) {
          cell.sub <- plot.meta[which(plot.meta[, sub.name] >= sub[1]), "cell"]
        } else {
          cell.sub <- plot.meta[which((plot.meta[, sub.name] >= sub[1]) & (plot.meta[, sub.name] < sub[2])), "cell"]
        }
      } else {
        cell.sub <- NULL
      }
      if (logical.connect == "and") {
        cell.left <- intersect(cell.left, cell.sub)
      } else if (logical.connect == "or") {
        cell.left <- union(cell.left, cell.sub)
      } else {
        stop(Sys.time(), " Unidentified logical.connect")
      }
    }
  }

  return(cell.left)
}




#'
#' constraintMatrix
#'
#' @name constraintMatrix
#'
#' @description
#' constraint FCS data by a provid cutoff
#'
#' @param x matrix
#' @param cutoff numeric. Cutoff of the constraint value
#' @param markers character. Markers used in the calculation of constraint model.
#' @param method character. the distance measure to be used.
#'    This must be one of "euclidean", "maximum", "manhattan",
#'    "canberra", "binary" or "minkowski".
#'
#' @export
#' @return  a matrix
#'
#' @examples
#'
#' mat <- matrix(runif(10000), nrow = 1000, ncol = 10)
#' colnames(mat) <- LETTERS[ seq_len(10)]
#' dim(mat)
#'
#' mat <- constraintMatrix(mat)
#' dim(mat)
#'
constraintMatrix <- function(x, cutoff = 0.99, markers = NULL, method = "euclidean") {

  if (!is.numeric(x)) stop(Sys.time(), " x must be a matrix ")

  if (is.null(markers)) markers <- colnames(x)
  if (!all(markers %in% colnames(x))) stop(Sys.time(), " markers must belong to the colnames of x ")

  sub <- abs(x[, markers])
  if (length(markers) > 1) {
    sub.mean <- colMeans(sub)
    d <- sapply(seq_len(nrow(sub)), function(aa) dist(rbind(sub[aa,], sub.mean), method = method))
  } else {
    sub.mean <- mean(sub)
    d <- sapply(seq_along(sub), function(aa) dist(rbind(sub[aa], sub.mean), method = method))
  }

  filter.mat <- x[order(d), ]
  filter.mat <- filter.mat[seq_len(floor(cutoff*dim(filter.mat)[1])), ]

  return(filter.mat)
}








