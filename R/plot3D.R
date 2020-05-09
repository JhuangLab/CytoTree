#'
#' Visualization of 3D data of FSPY
#'
#' @name plot3D
#'
#' @param object An FSPY object
#' @param item.use character. Items use to 3D plot, axes x and y and z must be numeric.
#' @param color.by character. Dot or mesh color by which character. It can be one of the column
#'     of plot.meta, or it can be just "density" (the default value).
#' @param order.by character. Order of color theme.
#' @param size numeric. size of the dot
#' @param angle numberic. angle of the plot
#' @param scale.y numeric. scale of y axis related to x- and z axis
#' @param category character. numeric or categorical
#' @param main character. title of the plot
#' @param color.theme vector. Color themes use in the plot.
#' @param ... options to pass on to the \code{\link[scatterplot3d]{scatterplot3d}} function.
#'
#' @import scatterplot3d
#' @importFrom grDevices rainbow colorRampPalette
#'
#' @export
#' @return gplots figure
#' @examples
#'
#' if (FALSE) {
#'
#'  plot3D(fspy, item.use = c("DC_2","DC_1","DC_3"), color.by = "stage",
#'         size = 0.5, angle = 60, color.theme = c("#FF99FF","#7A06A0","#FF3222"))
#'
#' }
#'
plot3D <- function(object,
                   item.use = c("PC1", "PC2", "PC3"),
                   color.by = "stage",
                   order.by = NULL,
                   size = 1,
                   angle = 60,
                   scale.y = 0.8,
                   category = "categorical",
                   main = "3D plot of FSPY",
                   color.theme = NULL,
                   ...) {

  # update and fetch plot meta information
  plot.meta <- fetchPlotMeta(object, verbose = FALSE)
  idx <- match(c(color.by, item.use), colnames(object@raw.data))
  idx <- idx[which(!is.na(idx))]
  if (length(idx) > 0) {
    sub <- as.data.frame(object@raw.data[which(object@meta.data$dowsample == 1), idx])
    colnames(sub) <- colnames(object@raw.data)[idx]
    plot.meta <- cbind(plot.meta, sub)
  }

  if ( !all(item.use %in% colnames(plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if ( !all(color.by %in% colnames(plot.meta)) ) stop(Sys.time(), " [ERROR] item.use is not in plot.meta of FSPY, please run updatePlotMeta first.")

  if (length(item.use) < 3) stop(Sys.time(), " [ERROR] item.use is less than two characters.")
  if (length(item.use) > 3) {
    warning(Sys.time(), " [WARNING] item.use is more than two characters. Only the first two will be used")
    item.use <- item.use[1:3]
  }
  if (length(color.by) > 1) {
    warning(Sys.time(), " [WARNING] color.by is more than one characters. Only the first one will be used")
    color.by <- color.by[1]
  }

  item.use.idx <- match(item.use, colnames(plot.meta))
  color.by.idx <- match(color.by, colnames(plot.meta))

  plot.data <- data.frame(plot.x = plot.meta[, item.use.idx[1]],
                          plot.y = plot.meta[, item.use.idx[2]],
                          plot.z = plot.meta[, item.use.idx[3]],
                          color.by = plot.meta[, color.by.idx])

  if ((length( unique(plot.data$color.by) ) > 256) & (category != "numeric")) {
    warning(Sys.time(), " [WARNING] color.by is categorical and has more than 50 elements. It will be used as numeric instead.")
    category = "numeric"
  }

  if (is.null(category)) {
    if (is.numeric(plot.data$color.by)) category="numeric" else category="categorical"
  }
  if (category == "categorical") {
    if (is.null(order.by)) {
      plot.data$color.by <- factor(plot.data$color.by)
    } else {
      plot.data$color.by <- factor(as.character(plot.data$color.by), levels = order.by)
    }
    if (is.null(color.theme)) {
      plot.data$color.by.3d <- factor(plot.data$color.by, labels = rainbow(length(levels(plot.data$color.by))))
    } else {
      plot.data$color.by.3d <- factor(plot.data$color.by, labels = colorRampPalette(color.theme)(length(levels(plot.data$color.by))) )
    }

  } else if (category == "numeric") {
    if (!is.numeric(plot.data$color.by)) plot.data$color.by <- as.numeric(factor(plot.data$color.by))

    if (is.null(color.theme)) {
      color.lib <- rainbow(102)
    } else {
      color.lib <- colorRampPalette(color.theme)(102)
    }
    plot.data$color.by.sd <- plot.data$color.by - min(plot.data$color.by)
    plot.data$color.by.3d <- color.lib[ ceiling( plot.data$color.by.sd/max(plot.data$color.by.sd) * 100 ) + 1 ]
  } else {
    warning(Sys.time(), " [WARNING] Unidentified parameters of category")
  }


  scatterplot3d(x = plot.data$plot.x, y = plot.data$plot.y, z = plot.data$plot.z,
                color = plot.data$color.by.3d,
                pch = 16, cex.symbols = size,
                scale.y = scale.y, angle = angle,
                xlab = item.use[1], ylab = item.use[2], zlab = item.use[3],
                main = main,
                col.axis = "#444444", col.grid = "#CCCCCC",
                ...)


}
