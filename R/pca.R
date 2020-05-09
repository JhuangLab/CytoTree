#'
#' Calculate principal components in FSPY
#'
#' @name runFastPCA
#'
#' @param object an FSPY object
#' @param center logical, a logical value indicating whether the variables
#'    should be shifted to be zero centered. Alternately, a vector
#'    of length equal the number of columns of x can be supplied.
#'    The value is passed to scale. See \code{\link[gmodels]{fast.prcomp}}
#' @param scale. logical, a logical value indicating whether the
#'    variables should be scaled to have unit variance before the
#'    analysis takes place. The default is FALSE for consistency
#'    with S, but in general scaling is advisable. Alternatively,
#'    a vector of length equal the number of columns of x can be supplied.
#'    The value is passed to scale. See \code{\link[gmodels]{fast.prcomp}}
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to \code{\link[gmodels]{fast.prcomp}} function
#'
#' @importFrom gmodels fast.prcomp
#' @seealso \code{\link[gmodels]{fast.prcomp}}
#'
#' @export
#' @return An FSPY object with PCA
#' @examples
#'
#' if (FALSE) {
#' fspy <- runFastPCA(fspy, verbose = TRUE)
#' }
#'
runFastPCA <- function(object, center = FALSE, scale. = TRUE,
                       verbose = FALSE, ...) {
  # PCA calculation
  if (verbose) message(Sys.time(), " [INFO] Calculating PCA.")
  if (length(which(object@meta.data$dowsample == 1)) < 10) stop(Sys.time, " [ERROR] Not enough cells, please run processingCluster and choose correct downsampleing.size paramter. ")
  mat <- object@log.data[which(object@meta.data$dowsample == 1), ]
  pca.obj <- fast.prcomp( t(mat), retx = TRUE, center = center, scale. = scale., ...)

  object@pca.sdev <- pca.obj$sdev
  object@pca.value <- pca.obj$rotation
  object@pca.scores <- pca.obj$x

  colnames(object@pca.value) <- paste0("PC_", 1:ncol(object@pca.value))

  if (verbose) message(Sys.time(), " [INFO] Calculating PCA completed. ")

  return(object)
}










