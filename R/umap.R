#'
#' Calculating UMAP
#'
#' @name runUMAP
#'
#' @description
#' Calculate Uniform Manifold Approximation and Projection in CYT
#'
#' @param object a CYT object
#' @param umap.config object of class umap.config. See \code{\link[umap]{umap}}.
#' @param n_neighbors numeric. Number of neighbors
#' @param dims numeric. Dim of umap, you can also change it in umap.config.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Options to pass on to the \code{\link[umap]{umap}} function
#'
#' @import umap
#' @seealso \code{\link[umap]{umap}}
#' @return A CYT object
#'
#' @export
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#'
#' cyt <- runUMAP(cyt, verbose = TRUE)
#' cyt <- runUMAP(cyt, n_neighbors = 20, verbose = TRUE)
#'
#' 
#'
#'
runUMAP <- function(object, umap.config = umap.defaults,
                    n_neighbors = 30, dims = 2, verbose = FALSE, ...) {
  if (verbose) message(Sys.time(), " Calculating Umap.")
  if (length(which(object@meta.data$dowsample == 1)) < 10) stop(Sys.time, " Not enough cells, please run processingCluster and choose correct downsampling.size paramter. ")
  mat <- as.matrix(object@log.data[which(object@meta.data$dowsample == 1), object@markers.idx])

  umap.config$n_neighbors <- n_neighbors
  umap.config$n_components <- dims
  umap.out <- umap(mat, config = umap.config, ...)
  object@umap.value <- umap.out$layout
  colnames(object@umap.value) <- paste0("UMAP_", seq_len(ncol(umap.out$layout)))
  rownames(object@umap.value) <- rownames(mat)

  if (verbose) message(Sys.time(), " Calculating Umap.")
  return(object)
}
