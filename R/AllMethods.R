#' CYT show method
#'
#' Prevents R from crashing by trying to print all slots of an CYT object  if a returned object is not stored in a variable.
#'
#' @param object An CYT object
#' @aliases show, show-method
#'
#' @export
#'
#' @docType methods
#' @method CytoTree show
#'
#' @keywords internal
#'
#' @return Cell number of CYT object
#'
#' @examples
#' if (FALSE) {
#'   cyt
#' }
#'
setMethod(
  f = "show",
  signature = "CYT",
  definition = function(object) {
    cat(
      "CYT Information:\n",
      "Input cell number:", nrow(object@raw.data), " cells \n",
      "Enroll marker number:", ncol(object@log.data), " markers \n",
      "Cells after downsampling:", sum(object@meta.data$dowsample), " cells \n"
    )
    invisible(NULL)
  }
)