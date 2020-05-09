#'
#' Apply gating on the matrix data
#'
#' @param x matrix
#' @param lower.gate vector. Gating parameter, the name of the vector is the marker name, and
#'    the value of the vector is the lower bound of gating cutoff.
#' @param upper.gate vector. Gating parameter, the name of the vector is the marker name, and
#'    the value of the vector is the upper bound of gating cutoff.
#'
#' @export
#'
#' @return a matrix
#'
#' @examples
#' par(mfrow=c(1,2))
#' x <- matrix(rnorm(200, 3, 1), nrow = 100, ncol = 2)
#' colnames(x) <- c("CD34", "CD43")
#' plot(x[, "CD34"], x[, "CD43"], main = "Before gating")
#'
#' lower.gate = c(CD34 = 2, CD43 = 3)
#' upper.gate = c(CD34 = 4, CD43 = 5)
#'
#' x <- gatingMatrix(x, lower.gate = lower.gate, upper.gate = upper.gate)
#' plot(x[, "CD34"], x[, "CD43"], main = "After gating")
#'
#' par(mfrow=c(1,1))
#'
#'
gatingMatrix <- function(x, lower.gate = NULL, upper.gate = NULL) {
  if (is.null(lower.gate) & is.null(upper.gate)) {
    warning(paste0(Sys.time(), " [WARNING] gating parameter is missing"))
    return(x)
  } else {
    if (!is.vector(lower.gate)) {
      warning(paste0(Sys.time(), " [WARNING] lower.gate parameter must be a vector"))
      return(x)
    } else if (!is.vector(upper.gate)) {
      warning(paste0(Sys.time(), " [WARNING] upper.gate parameter must be a vector"))
      return(x)
    } else {
      lname <- lower.gate[!names(lower.gate) %in% colnames(x)]
      uname <- upper.gate[!names(upper.gate) %in% colnames(x)]
      if (length(lname) > 0) {
        warning( paste0(Sys.time(), " [WARNING] some names in lower.gate is not in colnames of x") )
      }
      if (length(uname) > 0) {
        warning( paste0(Sys.time(), " [WARNING] some names in upper.gate is not in colnames of x") )
      }
      lname <- lower.gate[names(lower.gate) %in% colnames(x)]
      uname <- upper.gate[names(upper.gate) %in% colnames(x)]
      if (length(lname) == 0) {
        warning( paste0(Sys.time(), " [WARNING] names in lower.gate is not in colnames of x") )
        return(x)
      }
      if (length(uname) == 0) {
        warning( paste0(Sys.time(), " [WARNING] names in upper.gate is not in colnames of x") )
        return(x)
      }
      for (i in 1:length(lname)) {
        x <- x[which(x[, names(lname)[i]] > lname[i]), ]
        if (nrow(x) == 0) stop( paste0(Sys.time(), " [ERROR] lower.gate is out of bound") )
      }
      for (i in 1:length(uname)) {
        x <- x[which(x[, names(uname)[i]] < uname[i]), ]
        if (nrow(x) == 0) stop( paste0(Sys.time(), " [ERROR] upper.gate is out of bound") )
      }
      return(x)
    }
  }

}




