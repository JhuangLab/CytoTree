#'
#' subset FSPY object
#'
#' @name subsetFSPY
#'
#' @description This subsets an FSPY object by given a list of cells or cluster id.
#'     This function will subset all results without recalculating them, such as knn,
#'     PCA, tSNE, umap and pseudotime. For instance, you can choose recalculate PCA and
#'     tSNE and destiny scores by paramter recalculate.
#'
#' @param object An FSPY object
#' @param cells vector, Names of the cells to retain.
#' @param knn numeric. If is NA, the KNN will be equal to the knn number in the input FSPY object.
#' @param verbose logic. Whether to print calculation progress.
#'
#' @return An FSPY object
#'
#' @importFrom methods new
#'
#' @export
#'
#' @examples
#'
#' if (FALSE) {
#'
#' cells <- test.meta.data$cell[which(test.meta.data$stage == "D0")]
#' sub.fspy <- subsetFSPY(fspy, cells = cells)
#' sub.fspy
#'
#' }
#'
subsetFSPY <- function(object, cells = NULL,
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

  if (verbose) message(Sys.time(), " [INFO] Subset FSPY object.")
  object.new <- new("FSPY", raw.data = raw.data,
                    meta.data = meta.data,
                    log.data = log.data,
                    markers = object@markers, markers.idx = object@markers.idx)

  return(object.new)
}



