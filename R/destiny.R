#'
#' Calculate diffusion map in CYT
#'
#' @name runDiffusionMap
#'
#' @param object a CYT object
#' @param sigma.use numeric. Diffusion scale parameter of the Gaussian kernel. 
#'     One of '\code{local}',
#'     '\code{global}', a \code{\link[base]{numeric}} global sigma or a Sigmas object.
#'     When choosing '\code{global}', a global sigma will be calculated using find_sigmas
#'     (See \code{destiny}). A larger sigma might be necessary if the eigenvalues can not
#'    be found because of a singularity in the matrix. See \code{destiny}.
#' @param distance Distance measurement method applied to data or a distance matrix/dist.
#'    For the allowed values, see \code{destiny}
#' @param k numeric. By default is 30. \code{destiny} can be used to specify k.
#' @param density.norm logical. If TRUE, use density normalisation. See \code{destiny}
#' @param verbose logical. Whether to print calculation progress.
#' @param ... options to pass on to the \code{destiny}.
#'
#' @seealso \code{destiny}
#'
#' @import destiny
#'
#' @export
#' @return A CYT object
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' 
#' cyt <- runDiffusionMap(cyt, verbose = TRUE)
#' 
#'
runDiffusionMap <- function(object, sigma.use = NULL,
                            distance = c("euclidean", "cosine", "rankcor"),
                            k = 30,
                            density.norm = TRUE,  verbose = FALSE,
                            ...) {

  if (length(which(object@meta.data$dowsample == 1)) < 10) stop(Sys.time, " Not enough cells, please run processingCluster and choose correct downsampleing.size paramter. ")
  dm.data <- as.matrix(object@log.data[which(object@meta.data$dowsample == 1), object@markers.idx])

  if (verbose) message(Sys.time(), " Calculating Diffusion Map.")
  # Figure out sigma
  # this function refered to URD calcDM function.
  if (is.null(sigma.use)) {
    sigma.use <- find_sigmas(dm.data, verbose = FALSE)@optimal_sigma
    if (verbose) message(Sys.time(), " Destiny determined an optimal global sigma: ", round(sigma.use, digits=3))
  } else if (is.numeric(sigma.use)) {
    if (verbose) message(Sys.time(), " Using provided global sigma: ", round(sigma.use, digits=3))
  } else if (sigma.use == "local") {
    if (verbose) message(Sys.time(), " Using local sigma ")
  } else {
    sigma.use <- find_sigmas(dm.data, verbose = FALSE)@optimal_sigma
    warning(Sys.time(), " Invalid sigma value. Using an optimal global sigma instead.")
  }
  # Calculate the Diffusion Map
  distance <- match.arg(distance)
  dm.obj <- DiffusionMap(dm.data, sigma=sigma.use, k=k, density_norm = density.norm, distance=distance, ...)

  rownames(dm.obj@eigenvectors) <- rownames(dm.data)
  colnames(dm.obj@eigenvectors) <- paste0("DC_", seq_len(ncol(dm.obj@eigenvectors)))
  rownames(dm.obj@transitions) <- rownames(dm.data)
  colnames(dm.obj@transitions) <- rownames(dm.data)


  # Load diffusion map into the Dropseq object
  object@dm <- dm.obj

  if (verbose) message(Sys.time(), " Calculating Diffusion Map completed")

  return(object)
}





