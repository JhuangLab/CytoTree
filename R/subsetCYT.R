#'
#' subset CYT object
#'
#' @name subsetCYT
#'
#' @description This subsets an CYT object by given a list of cells or cluster id.
#'     This function will subset all results without recalculating them, such as knn,
#'     PCA, tSNE, umap and pseudotime. For instance, you can choose recalculate PCA and
#'     tSNE and destiny scores by paramter recalculate.
#'
#' @param object An CYT object
#' @param cells vector, Names of the cells to retain.
#' @param knn numeric. If is NA, the KNN will be equal to the knn number in the input CYT object.
#' @param verbose logic. Whether to print calculation progress.
#'
#' @return An CYT object
#'
#' @importFrom methods new
#'
#' @export
#'
#' @examples
#'
#' cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
#' cyt <- readRDS(file = cyt.file)
#'
#' meta.data <- fetchPlotMeta(cyt)
#' cells <- meta.data$cell[which(meta.data$stage == "D0")]
#' sub.cyt <- subsetCYT(cyt, cells = cells)
#' sub.cyt
#'
#' 
#'
subsetCYT <- function(object, cells = NULL,
                       knn = NA,
                       verbose = FALSE) {
  if (is.null(cells)) {
    warning(Sys.time(), " [WARNING] cells must be provided.")
    cells <- rownames(object@raw.data)
  }
  # Make sure all cells are actually in the object
  cells.keep <- intersect(cells, rownames(object@raw.data))

  raw.data <- object@raw.data[cells.keep, ]
  log.data <- object@log.data[cells.keep, ]
  meta.data <- object@meta.data[cells.keep, ]

  if (verbose) message(Sys.time(), " [INFO] Subset CYT object.")
  object.new <- new("CYT", raw.data = raw.data,
                    meta.data = meta.data,
                    log.data = log.data,
                    markers = object@markers, markers.idx = object@markers.idx)

  return(object.new)
}



