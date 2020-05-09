
#'
#' plot MST of FSPY
#'
#' @name plotTree
#'
#' @param object an FSPY object
#' @param cex.size numeric. size cex of the dot
#' @param color.by numeric. size color theme of the dot
#' @param size.by numeric. size theme of the dot
#' @param as.tree logical. Whether to show node as tree
#' @param root.id numeric. Root id of the tree, if as.tree is TRUE
#' @param show.node.name logical. whether to show node name
#'
#' @export
#' @importFrom igraph layout_as_tree layout.kamada.kawai as_data_frame
#' @return ggplot2 figure
#'
#' @examples
#'
#' if (FALSE) {
#'
#' plotTree(fspy)
#'
#' plotTree(fspy, show.node.name = T)
#'
#' plotTree(fspy, color.by = "CD43", show.node.name = T, cex.size = 1) +
#'     scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222"))
#'
#' plotTree(fspy, color.by = "D0.percent", show.node.name = T, cex.size = 1) +
#'     scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222"))
#'
#' plotTree(fspy, color.by = "D2.percent", show.node.name = T, cex.size = 1) +
#'     scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222"))
#'
#' plotTree(fspy, color.by = "pseudotime", cex.size = 1) +
#'     scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))
#'
#' }
#'
plotTree <- function(object,
                     cex.size = 1,
                     color.by = "cell.number",
                     size.by = "cell.number",
                     as.tree = FALSE,
                     root.id = NULL,
                     show.node.name = FALSE) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  if (is.null(object@network)) stop(Sys.time(), " [ERROR] network is missing, please run runCluster first!")

  mst <- object@network$mst

  # update plot meta information
  node.attr <- fetchClustMeta(object, verbose = FALSE)

  edge.attr <- igraph::as_data_frame(mst)

  ##### layout
  if (as.tree) {
    if (is.null(root.id)) {
      root.id = node.attr$cluster[which(node.attr$pseudotime == min(node.attr$pseudotime))]
    }
    l <- igraph::layout_as_tree(mst, root = root.id)
  } else {
    l <- igraph::layout.kamada.kawai(mst)
  }
  colnames(l) <- c("pos.x", "pos.y")
  node.attr <- cbind(node.attr, l)

  size.by.idx <- match(size.by, colnames(node.attr))
  color.by.idx <- match(color.by, colnames(node.attr))

  edge.attr$from.x <- node.attr$pos.x[match(edge.attr$from, node.attr$cluster)]
  edge.attr$from.y <- node.attr$pos.y[match(edge.attr$from, node.attr$cluster)]
  edge.attr$to.x <- node.attr$pos.x[match(edge.attr$to, node.attr$cluster)]
  edge.attr$to.y <- node.attr$pos.y[match(edge.attr$to, node.attr$cluster)]

  color.tree <- node.attr[, color.by.idx]
  size.tree <- node.attr[, size.by.idx]

  gg <- ggplot()
  gg <- gg + geom_segment(mapping = aes(x = edge.attr$from.x, y = edge.attr$from.y, xend = edge.attr$to.x, yend = edge.attr$to.y))
  gg <- gg + geom_point(mapping = aes(x = node.attr$pos.x, y = node.attr$pos.y, color = color.tree, size = size.tree))
  gg <- gg + scale_size(range = c(0, 6) * cex.size)
  gg <- gg + labs(color = color.by)
  gg <- gg + labs(size = size.by)

  if (show.node.name) gg <- gg + geom_text(aes(x = node.attr$pos.x, y = node.attr$pos.y, label = node.attr$cluster ), check_overlap = TRUE, size = 3 * cex.size)
  gg <- gg + theme_void()
  gg <- gg + labs(x = "", y = "", title = paste0("Tree plot, color.by: ", color.by, ", size.by: ", size.by))

  return(gg)

}


#'
#' plot MST pie of FSPY
#'
#' @name plotPieTree
#'
#' @param object an FSPY object
#' @param cex.size numeric. size cex of the dot
#' @param size.by.cell.number logical. Whether to size node by cell number
#' @param as.tree logical. Whether to show node as tree
#' @param root.id numeric. Root id of the tree, if as.tree is TRUE
#' @param show.node.name logical. whether to show node name
#'
#' @export
#' @importFrom igraph layout_as_tree layout.kamada.kawai as_data_frame
#' @import scatterpie
#' @return ggplot2 figure
#'
#' @examples
#'
#' if (FALSE) {
#'
#' # Runs only have two or more stages
#' plotPieTree(fspy, cex.size = 1, size.by.cell.number = T) +
#'    scale_fill_manual(values = c("#00599F","#FF3222","#009900",
#'                                 "#FF9933","#FF99FF","#7A06A0"))
#' }
#'
plotPieTree <- function(object,
                        cex.size = 2,
                        size.by.cell.number = TRUE,
                        as.tree = FALSE,
                        root.id = NULL,
                        show.node.name = FALSE) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  if (is.null(object@network)) stop(Sys.time(), " [ERROR] network is missing, please run runCluster first!")
  if (length(unique(object@meta.data$stage)) <= 1) stop(Sys.time(), " [ERROR] plotPieTree only fits elements in stage over 2!")

  mst <- object@network$mst

  # update plot meta information
  node.attr <- fetchClustMeta(object, verbose = FALSE)

  edge.attr <- igraph::as_data_frame(mst)
  pos.x <- pos.y <- cluster <- cell.number.percent <- NULL
  ##### layout
  if (as.tree) {
    if (is.null(root.id)) {
      root.id = node.attr$cluster[which(node.attr$pseudotime == min(node.attr$pseudotime))]
    }
    l <- igraph::layout_as_tree(mst, root = root.id)
  } else {
    l <- igraph::layout.kamada.kawai(mst)
  }
  colnames(l) <- c("pos.x", "pos.y")
  node.attr <- cbind(node.attr, l)

  plot.cols <- paste0(unique(object@meta.data$stage), ".percent")

  edge.attr$from.x <- node.attr$pos.x[match(edge.attr$from, node.attr$cluster)]
  edge.attr$from.y <- node.attr$pos.y[match(edge.attr$from, node.attr$cluster)]
  edge.attr$to.x <- node.attr$pos.x[match(edge.attr$to, node.attr$cluster)]
  edge.attr$to.y <- node.attr$pos.y[match(edge.attr$to, node.attr$cluster)]

  gg <- ggplot()
  gg <- gg + geom_segment(mapping = aes(x = edge.attr$from.x, y = edge.attr$from.y, xend = edge.attr$to.x, yend = edge.attr$to.y))
  if (size.by.cell.number) {
    gg <- gg + geom_scatterpie(aes(x = pos.x, y = pos.y, group = cluster, r = cell.number.percent*cex.size),
                               data = node.attr, cols = plot.cols, color=NA) + coord_equal()
  } else {
    gg <- gg + geom_scatterpie(aes(x = pos.x, y = pos.y, group = cluster, r = 0.1*cex.size),
                               data = node.attr, cols = plot.cols, color=NA) + coord_equal()
  }

  gg <- gg + theme_void()
  gg <- gg + labs(x = "", y = "", title = paste0("Pie Tree Plot, size by cell.number"))


  return(gg)

}
















