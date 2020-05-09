#'
#' plot Pseudotime density of FSPY
#'
#' @name plotPseudotimeDensity
#'
#' @param object an FSPY object
#' @param color.by character.
#' @param main character. Title of the plot
#' @param adjust numeric. A multiplicate bandwidth adjustment.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @export
#' @return ggplot2 figure
#'
#' @examples
#'
#' if (FALSE) {
#'
#' plotPseudotimeDensity(fspy)
#'
#' plotPseudotimeDensity(fspy, adjust = 1)
#' plotPseudotimeDensity(fspy, adjust = 2)
#'
#' plotPseudotimeDensity(fspy, adjust = 2) +
#'   scale_color_manual(values = c("#00599F","#FF3222","#009900",
#'                                 "#FF9933","#FF99FF","#7A06A0"))
#'
#' }
#'
plotPseudotimeDensity <- function(object, color.by = "stage",
                                  main = "Density of pseudotime",
                                  adjust = 0.5,
                                  plot.theme = theme_bw()) {
  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  object <- updatePlotMeta(object, verbose = FALSE)

  pseudotime = Pseudotime = Signal = Marker = Stage = NULL

  # checking items
  if ( !all("pseudotime" %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] pseudotime is not in plot.meta of FSPY, please run Pseudotime first.")

  if ( !all(color.by %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  item.use.idx <- match("pseudotime", colnames(object@plot.meta))
  color.by.idx <- match(color.by, colnames(object@plot.meta))

  plot.data <- data.frame(pseudotime = object@plot.meta[, item.use.idx[1]],
                          color.by = object@plot.meta[, color.by.idx])

  if (length(unique(plot.data$color.by)) > 50) {
    message(Sys.time(), " [INFO] color.by is a numeric vector and has over 50 elements")
    plot.data$color.by <- 1
  }
  if (!is.factor(plot.data$color.by)) plot.data$color.by <- as.factor(as.character(plot.data$color.by))

  gg <- ggplot(plot.data, aes(x=pseudotime, colour = color.by))
  gg <- gg + geom_density(adjust = adjust)
  gg <- gg + plot.theme
  gg <- gg + labs(color = color.by)

  gg <- gg + labs(title = paste0(main))

  return(gg)
}


#'
#' plotPseudotimeTraj
#'
#' @name plotPseudotimeTraj
#'
#' @param object An FSPY object
#' @param markers character. Markers used in the calculation progress
#' @param cutoff numeric. Cutoff of trajectory value
#' @param size numeric. Size of the dot
#' @param alpha numeric. Transparency (0-1) of the dot, default is 1.
#' @param print.curve logical. Whether to perform curve fitting
#' @param var.cols logical. Whether to plot stage
#' @param plot.theme themes from \code{ggplot2}
#'
#' @importFrom stats predict
#'
#' @export
#' @return ggplot2 figure
#'
#' @examples
#'
#' if (FALSE) {
#'
#' plotPseudotimeTraj(fspy)
#' plotPseudotimeTraj(fspy, print.curve = F)
#' plotPseudotimeTraj(fspy, var.cols = T)
#'
#' plotPseudotimeTraj(fspy) +
#'    scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))
#'
#' plotPseudotimeTraj(fspy, markers = c("CD43", "CD34")) +
#' scale_colour_gradientn(colors = c("#F4D31D", "#FF3222","#7A06A0"))
#'
#' }
#'
plotPseudotimeTraj <- function(object,
                               cutoff = -1,
                               markers = NULL,
                               size = 0.5,
                               alpha = 0.6,
                               print.curve = TRUE,
                               var.cols = FALSE,
                               plot.theme = theme_bw()) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  object <- updatePlotMeta(object, verbose = FALSE)

  # checking items
  if ( !all("pseudotime" %in% colnames(object@plot.meta)) ) stop(Sys.time(), " [ERROR] pseudotime is not in plot.meta of FSPY, please run Pseudotime first.")

  # checking items
  if (cutoff > 0) {
    if ( !all(c("traj.value","traj.value.log") %in% colnames(object@plot.meta)) ) {
      message(Sys.time(), " [INFO] traj.value is not in plot.meta of FSPY, please run runWalk first.")
    }
  }
  if (is.null(markers)) markers <- object@markers
  if (length(markers) > 30) {
    warning(Sys.time(), " [WARNING] only the first 30 markers will be plot")
    markers <- markers[1:30]
  }

  plot.data <- NULL
  plot.meta <- object@plot.meta
  Pseudotime <- Signal <- Marker <- Stage <- NULL
  for (i in 1:length(markers)) {
    sub <- data.frame(Pseudotime = plot.meta$pseudotime,
                      IsRoot = plot.meta$is.root.cells,
                      IsLeaf = plot.meta$is.leaf.cells,
                      Marker = markers[i],
                      Signal = object@raw.data[which(object@meta.data$dowsample == 1), markers[i]],
                      Stage = plot.meta$stage,
                      TrajValue = plot.meta$traj.value,
                      LogTrajValue = plot.meta$traj.value.log)
    idx <- which( (sub$LogTrajValue > cutoff)  )
    if (sum(sub$IsRoot == 1) > 0) {
      idx.2 <- which(sub$IsRoot == 1)
    } else {
      idx.2 <- NULL
    }
    idx <- union(idx, idx.2)
    if (sum(sub$IsLeaf == 1) > 0) {
      idx.3 <- which(sub$IsLeaf == 1)
    } else {
      idx.3 <- NULL
    }
    idx <- union(idx, idx.3)
    sub <- sub[idx, ]
    plot.data <- rbind(plot.data, sub)
  }

  gg <- ggplot(plot.data, aes(x=Pseudotime, y=Signal, color = Pseudotime)) + geom_point(size = size, alpha = alpha)
  gg <- gg + plot.theme

  if (var.cols) {
    gg <- gg + facet_grid(rows = vars(Marker), cols = vars(Stage))
  } else {
    gg <- gg + facet_grid(rows = vars(Marker))
  }

  if (print.curve) {
    gg <- gg + geom_smooth(color="black", method="loess", se = FALSE)
  }

  return(gg)
}


#'
#' plotMarkerDensity
#'
#' @name plotMarkerDensity
#'
#' @param object An FSPY object
#' @param markers character. Markers used in the calculation progress
#' @param cutoff numeric. Cutoff of trajectory value
#' @param adjust numeric. Transparency (0-1) of the dot, default is 1.
#' @param plot.theme themes from \code{ggplot2}
#'
#' @importFrom stats predict
#'
#' @export
#' @return ggplot2 figure
#'
#' @examples
#'
#' if (FALSE) {
#'
#' plotMarkerDensity(fspy)
#' plotMarkerDensity(fspy, adjust = 1)
#'
#' }
#'
plotMarkerDensity <- function(object,
                              cutoff = -1,
                              markers = NULL,
                              adjust = 0.5,
                              plot.theme = theme_bw()) {

  if (missing(object)) stop(Sys.time(), " [ERROR] object is missing")
  object <- updatePlotMeta(object, verbose = FALSE)

  # checking items
  if (cutoff > 0) {
    if ( !all(c("traj.value","traj.value.log") %in% colnames(object@plot.meta)) ) {
      message(Sys.time(), " [INFO] traj.value is not in plot.meta of FSPY, please run runWalk first.")
    }
  }
  if (is.null(markers)) markers <- object@markers
  if (length(markers) > 30) {
    warning(Sys.time(), " [WARNING] only the first 30 markers will be plot")
    markers <- markers[1:30]
  }

  plot.data <- NULL
  Pseudotime = IsRoot = IsLeaf = Marker = Signal = Stage = TrajValue = LogTrajValue = NULL
  plot.meta <- fetchPlotMeta(object, verbose = FALSE)
  for (i in 1:length(markers)) {
    sub <- data.frame(Pseudotime = plot.meta$pseudotime,
                      IsRoot = plot.meta$is.root.cells,
                      IsLeaf = plot.meta$is.leaf.cells,
                      Marker = markers[i],
                      Signal = object@raw.data[which(object@meta.data$dowsample == 1), markers[i]],
                      Stage = plot.meta$stage,
                      TrajValue = plot.meta$traj.value,
                      LogTrajValue = plot.meta$traj.value.log)
    idx <- which( (sub$LogTrajValue > cutoff)  )
    if (sum(sub$IsRoot == 1) > 0) {
      idx.2 <- which(sub$IsRoot == 1)
    } else {
      idx.2 <- NULL
    }
    idx <- union(idx, idx.2)
    if (sum(sub$IsLeaf == 1) > 0) {
      idx.3 <- which(sub$IsLeaf == 1)
    } else {
      idx.3 <- NULL
    }
    idx <- union(idx, idx.3)
    sub <- sub[idx, ]
    plot.data <- rbind(plot.data, sub)
  }

  gg <- ggplot(plot.data, aes(x=Signal, color = Marker)) + geom_density(adjust = adjust)
  gg <- gg + plot.theme
  gg <- gg + facet_grid(rows = vars(Stage), cols = vars(Marker))
  return(gg)
}




