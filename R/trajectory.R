
#'
#' Walk between root cells and leaf cells
#'
#' @name runWalk
#'
#' @description Walk between root cells and leaf cells
#'
#' @param object An FSPY object
#' @param mode character. Specifies how igraph should interpret the supplied matrix.
#'    Possible values are: undirected, directed, upper, lower, max, min, plus. By 
#'    default is undirected.
#' @param max.run.forward numeric. Maximum cycles of forward walk.
#' @param backward.walk logical. Whether to run backward walk.
#' @param max.run.backward numeric. Maximum cycles of backward walk.
#' @param verbose logical. Whether to print calculation progress.
#' @param ... Parameters passing to calculation function.
#'
#' @importFrom igraph graph.adjacency simplify shortest_paths
#'
#' @return An FSPY object
#'
#' @export
#'
#' @examples
#'
#' if (FALSE) {
#'   fspy <- runWalk(fspy, verbose = TRUE)
#'   fspy <- runWalk(fspy, backward.walk = FALSE, verbose = TRUE)
#' }
#'
#'
#'
runWalk <- function(object, mode = c("undirected", "directed", "max", "min", "upper", "lower", "plus"),
                    max.run.forward = 20,
                    backward.walk = FALSE, max.run.backward = 20,
                    verbose = FALSE, ...) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing.")

  if (dim(object@knn.index)[1] == 0) stop(Sys.time(), " [ERROR] KNN information is missing in FSPY, please run runKNN first.")

  if (verbose) message(Sys.time(), " [INFO] Calculating walk between root.cells and leaf.cells .")

  if (!"pseudotime" %in% colnames(object@meta.data)) stop(Sys.time(), " [INFO] Pseudotime exists in meta.data, it will be replaced.")

  knn.index <- object@knn.index

  adj <- matrix(0, nrow(knn.index), nrow(knn.index))
  rownames(adj) <- colnames(adj) <- rownames(knn.index)
  pseudotime <- object@meta.data$pseudotime[which(object@meta.data$seed.pseudotime == 1)]

  # generating a adjacency matrix by nearest neighbors
  if (verbose) message(Sys.time(), " [INFO] Generating an adjacency matrix.")
  for(i in seq_len(nrow(knn.index))) {
    idx <- knn.index[i,][ pseudotime[knn.index[i,]] > pseudotime[i]  ]
    adj[i, rownames(knn.index)[idx]] <- 1
  }
  mode <- match.arg(mode)
  g <- igraph::graph.adjacency(adj, mode = mode, ...)
  # remove self loops
  g <- simplify(g)

  if (verbose) message(Sys.time(), " [INFO] Walk forward.")
  root.cells <- object@root.cells[object@root.cells %in% rownames(knn.index)]
  leaf.cells <- object@leaf.cells[object@leaf.cells %in% rownames(knn.index)]
  if (length(root.cells) >= max.run.forward ) {
    root.cells <- as.character(sample(root.cells, max.run.forward))
  } else {
    warning(Sys.time(), " [WARNING] max.run.forward is too large.")
  }
  # run forward
  walk.forward <- suppressWarnings(lapply(as.character(root.cells), function(x) shortest_paths(g, from = x, to = as.character(leaf.cells))$vpath ))

  # run run backward
  if (backward.walk) {
    if (verbose) message(Sys.time(), " [INFO] Walk backward.")
    leaf.cells <- object@leaf.cells[object@leaf.cells %in% rownames(knn.index)]
    if (length(leaf.cells) >= max.run.backward ) {
      leaf.cells <- as.character(sample(leaf.cells, max.run.backward))
    } else {
      warning(Sys.time(), " [WARNING] max.run.backward is too large.")
    }
    walk.backward <- suppressWarnings(lapply(leaf.cells, function(x) shortest_paths(g, from = x, to = object@root.cells)$vpath ))
  } else {
    walk.backward <- NULL
    max.run.backward <- 0
  }

  object@meta.data$traj.value <- 0
  object@meta.data$traj.value.log <- 0

  cell.info <- unlist(c(walk.forward, walk.backward))
  cell.info <- as.data.frame(table(names(cell.info)))

  object@meta.data$traj.value[match(cell.info$Var1, object@meta.data$cell)] <- cell.info$Freq / (max.run.forward + max.run.backward)
  object@meta.data$traj.value[match(object@root.cells, object@meta.data$cell)] = 0
  object@meta.data$traj.value[match(object@leaf.cells, object@meta.data$cell)] = 0

  object@walk <- list(max.run.forward = max.run.forward,
                      max.run.backward = max.run.backward)

  object@meta.data$traj.value.log <- log10(object@meta.data$traj.value + 1)

  if (verbose) message(Sys.time(), " [INFO] Calculating walk completed.")

  return(object)
}





