#'
#' Update plot meta information of FSPY
#'
#' @name updatePlotMeta
#'
#' @param object An FSPY object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#' @return An FSPY object
#'
#' @examples
#'
#' if (FALSE) {
#' fspy <- updatePlotMeta(fspy)
#' }
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

  if (verbose) message(Sys.time(), " [INFO] Columns can be used in plot2D and plot3D: ", paste(colnames(plot.meta), collapse = " "))

  return(object)
}

#'
#' Update clusters' meta information of FSPY
#'
#' @name updateClustMeta
#'
#' @param object An FSPY object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @export
#' @return An FSPY object
#'
#' @importFrom stats aggregate
#'
#' @examples
#'
#' if (FALSE) {
#' fspy <- updateClustMeta(fspy)
#' }
#'
updateClustMeta <- function(object, verbose = TRUE) {

  # Generating tree meta information
  plot.data <- fetchPlotMeta(object, verbose = FALSE)
  plot.data <- cbind(plot.data, object@raw.data[which(object@meta.data$dowsample == 1), ])

  if (length(unique(plot.data$stage)) > 1) {
    cell.count <- table(plot.data[, match(c("cluster.id", "stage"), colnames(plot.data)) ])
    cell.count<- t(sapply(1:dim(cell.count)[1], function(x) cell.count[x, ]))
    cell.total.number <- rowSums(cell.count)
    cell.total.number.percent <- cell.total.number/sum(cell.total.number)
    cell.percent<- t(sapply(1:dim(cell.count)[1], function(x) cell.count[x,]/sum(cell.count[x,]) ))
    colnames(cell.percent) <- paste0(colnames(cell.count), ".percent")
  } else {
    cell.count <- table(plot.data[, match(c("cluster.id", "stage"), colnames(plot.data)) ])
    cell.count<- t(sapply(1:dim(cell.count)[1], function(x) cell.count[x, ]))
    cell.count <- as.data.frame(t(cell.count))
    cell.total.number <- rowSums(cell.count)
    cell.total.number.percent <- cell.total.number/sum(cell.total.number)
    cell.percent<- t(sapply(1:dim(cell.count)[1], function(x) cell.count[x,]/sum(cell.count[x,]) ))
    cell.percent <- as.data.frame(t(cell.percent))
    colnames(cell.count) <- paste0(unique(plot.data$stage))
    colnames(cell.percent) <- paste0(unique(plot.data$stage), ".percent")
  }

  idx.redim <- match(c(colnames(object@raw.data), "pseudotime", "traj.value", "traj.value.log"), colnames(plot.data))
  idx.redim <- unique(idx.redim)
  tree.meta <- stats::aggregate(plot.data[, idx.redim], list(cluster = plot.data[, "cluster.id"]), mean)
  tree.meta.1 <- data.frame(cell.count,
                           cell.number = cell.total.number,
                           cell.number.percent = cell.total.number.percent,
                           cell.percent)
  tree.meta <- cbind(tree.meta, tree.meta.1)
  tree.meta$branch.id <- plot.data$branch.id[match(tree.meta$cluster, plot.data$cluster.id)]

  object@tree.meta <- tree.meta
  if (verbose) message(Sys.time(), " [INFO] Columns can be used in plotTree: ", paste(colnames(tree.meta), collapse = " "))

  return(object)
}


#'
#' Fetching plot metadata of FSPY
#'
#' @name fetchPlotMeta
#'
#' @param object An FSPY object
#' @param markers vector. Makers fetched from expression matrix
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return a data.frame containing meta information for visualization
#'
#' @export
#'
#' @examples
#'
#' if (FALSE) {
#' plot.data <- fetchPlotMeta(fspy)
#' head(plot.data)
#'
#' plot.data <- fetchPlotMeta(fspy, markers = c("CD43", "CD34"))
#' head(plot.data)
#' }
#'
#'
#'
fetchPlotMeta <- function(object, markers = NULL, verbose = FALSE) {

    # update and fetch plot meta information
  object <- updatePlotMeta(object, verbose = FALSE)
  plot.meta <- object@plot.meta
  idx <- match(markers, colnames(object@raw.data))
  idx <- idx[which(!is.na(idx))]
  if (length(idx) > 0) {
    plot.meta <- cbind(plot.meta, object@raw.data[which(object@meta.data$dowsample == 1), idx])
  }


  return(plot.meta)
}

#'
#' Fetching clusters' metadata of FSPY
#'
#' @name fetchClustMeta
#'
#' @param object An FSPY object
#' @param verbose logical. Whether to print calculation progress.
#'
#' @return a data.frame containing clustering information for visualization
#'
#' @export
#'
#' @examples
#'
#' if (FALSE) {
#' clust.data <- fetchClustMeta(fspy)
#' head(clust.data)
#' }
#'
#'
fetchClustMeta <- function(object, verbose = FALSE) {

  object <- updateClustMeta(object, verbose = verbose)

  return(object@tree.meta)
}


#'
#' Fetching cellls of FSPY
#'
#' @name fetchCell
#'
#' @param object An FSPY object
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
#' if (FALSE) {
#' cell.fetch <- fetchCell(fspy, traj.value.log = 0.01)
#' cell.fetch <- fetchCell(fspy, stage = c("D0", "D10"))
#' cell.fetch <- fetchCell(fspy, stage = c("D0", "D10"), traj.value.log = 0.01,
#'                         logical.connect = "or")
#' }
#'
fetchCell <- function(object, logical.connect = "or", verbose = FALSE, ... ) {

  object <- updatePlotMeta(object, verbose = verbose)

  param.list <- list(...)
  plot.meta <- object@plot.meta

  cell.left <- NULL

  if (length(param.list) > 0) {
    param.list <- param.list[names(param.list) %in% colnames(plot.meta)]
    for (i in 1:length(param.list)) {
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
        stop(Sys.time(), " [ERROR] Unidentified logical.connect")
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
#' colnames(mat) <- LETTERS[1:10]
#' dim(mat)
#'
#' mat <- constraintMatrix(mat)
#' dim(mat)
#'
constraintMatrix <- function(x, cutoff = 0.99, markers = NULL, method = "euclidean") {

  if (!is.numeric(x)) stop(Sys.time(), " [ERROR] x must be a matrix ")

  if (is.null(markers)) markers <- colnames(x)
  if (!all(markers %in% colnames(x))) stop(Sys.time(), " [ERROR] markers must belong to the colnames of x ")

  sub <- abs(x[, markers])
  if (length(markers) > 1) {
    sub.mean <- colMeans(sub)
    d <- sapply(1:nrow(sub), function(aa) dist(rbind(sub[aa,], sub.mean), method = method))
  } else {
    sub.mean <- mean(sub)
    d <- sapply(1:length(sub), function(aa) dist(rbind(sub[aa], sub.mean), method = method))
  }

  filter.mat <- x[order(d), ]
  filter.mat <- filter.mat[1:floor(cutoff*dim(filter.mat)[1]), ]

  return(filter.mat)
}








