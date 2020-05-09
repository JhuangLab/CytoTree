#'
#' Calculate k-nearest neighbors of FSPY
#'
#' @name runKNN
#'
#' @description Calculates and stores a k-nearest neighbor graph based on Euclidean
#'    distance with (KMKNN) algorithm using log-transformed signaling matrix of
#'    flow cytometry data. The base function are base on \code{\link[BiocNeighbors]{findKNN}}.
#'
#' @param object an FSPY object
#' @param given.mat matrix. Given matrix to run knn
#' @param knn numeric. Number of k-nearest neighbors.
#' @param knn.replace logic. Whether to replace knn in FSPY object
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[BiocNeighbors]{findKNN}} function
#'
#' @seealso \code{\link[BiocNeighbors]{findKNN}}
#'
#' @return An FSPY object with knn, knn.index and knn.distance information.
#'
#' @import BiocNeighbors
#'
#' @export
#'
#' @examples
#' if (FALSE) {
#'
#' fspy <- runKNN(fspy)
#'
#' }
#'
#'
runKNN <- function(object,
                   given.mat = NULL,
                   knn = 30,
                   knn.replace = TRUE, 
                   verbose = FALSE, ...) {

  if (isTRUE(object@knn > 0) & !(knn.replace)) {
    if (verbose) message(Sys.time(), " [INFO] Using knn in FSPY object: ", object@knn )
  } else if ( isTRUE(object@knn > 0) & (knn.replace) ) {
    if (verbose) message(Sys.time(), " [INFO] Using knn provided in this function: ", knn )
    object@knn <- knn
  } else {
    object@knn <- knn
  }

  if (length(which(object@meta.data$dowsample == 1)) < 10) {
    stop(Sys.time, " [ERROR] Not enough cells, please run processingCluster and choose correct downsampleing.size paramter. ")
  }

  if (is.null(given.mat)) {
    mat <- object@log.data[which(object@meta.data$seed.pseudotime == 1), ]
  } else {
    if (nrow(given.mat) != nrow(object@log.data[which(object@meta.data$seed.pseudotime == 1), ])) {
      stop(Sys.time, " [ERROR] Invalid given.mat ")
    } else {
      mat <- given.mat
    }
  }

  if (verbose) message(paste0(Sys.time(), " [INFO] Calculating KNN " ) )
  fout <- suppressWarnings(findKNN(mat, k = object@knn, ...))

  rownames(fout$index) <- rownames(mat)
  rownames(fout$distance) <- rownames(mat)

  object@knn = knn
  object@knn.index = fout$index
  object@knn.distance = fout$distance

  if (verbose) message(Sys.time(), " [INFO] Calculating KNN completed. ")
  return(object)
}

