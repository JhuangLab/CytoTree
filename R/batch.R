#'
#' correctBatchCYT
#'
#' @description
#' Remove batch effect in CYT object
#' @param object a CYT object
#' @param batch vector. Batch covariate (only one batch allowed)
#' @param par.prior logical. TRUE indicates parametric adjustments will be used,
#'    FALSE indicates non-parametric adjustments will be used.
#' @param mean.only logical. FALSE If TRUE ComBat only corrects the mean of the batch
#'    effect (no scale adjustment)
#' @param verbose logical. Whether to show log information
#' @param ... Parameters passing to \code{\link[sva]{ComBat}} function
#'
#' @seealso \code{\link[BiocNeighbors]{findKNN}}
#'
#' @importFrom sva ComBat
#' @export
#'
#' @return A CYT object with corrected batch effects
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#' plot.meta <- fetchPlotMeta(cyt)
#' batch <- as.numeric(plot.meta$stage)
#' cyt <- correctBatchCYT(cyt, batch = batch)
#' 
#'
correctBatchCYT <- function(object, batch = NULL, par.prior = TRUE,
                            mean.only = TRUE, verbose = FALSE, ...) {
  log.data <- object@log.data[, object@markers.idx]

  # correct batch effect using ComBat
  suppressMessages(
    log.data.combat <- sva::ComBat(dat = t(log.data), batch = batch,
                                   par.prior = par.prior,
                                   mean.only = mean.only, ...))

  object@log.data[, object@markers.idx] <- t(log.data.combat)
  return(object)
}

