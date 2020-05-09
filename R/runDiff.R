#'
#' Calculate differential expression markers
#'
#' @name runDiff
#'
#' @description Calculating differentially expressed markers
#'
#' @param object an FSPY object
#' @param branch.id vector. Branch ids use to run differentially expressed markers
#' @param branch.id.2 vector. Branch ids use to run differentially expressed
#'    markers in compare with branch.id
#' @param verbose logic. Whether to print calculation progress.
#'
#' @seealso  \code{bulidTree}
#'
#' @return An FSPY object with cluster.id in meta.data
#'
#' @import limma
#' @importFrom stringr str_replace_all fixed
#' @importFrom stats model.matrix
#'
#' @export
#' @return a data.frame with differential expressed markers
#'
#' @examples
#'
#' if (FALSE) {
#'
#' DEG.table <- runDiff(fspy)
#'
#' }
#'
#'
#'
runDiff <- function(object, branch.id = NULL, branch.id.2 = NULL, verbose = FALSE) {

  if (verbose) message(Sys.time(), " [INFO] Calculating differentially expressed markers.")
  if (missing(object)) stop(Sys.time(), " [ERROR] FSPY object is missing.")
  if (!"branch.id" %in% colnames(object@meta.data)) stop(Sys.time(), " [ERROR] branch.id is missing, please run buildTree first.")

  all.branch.ids <- unique(object@meta.data$branch.id)

  total.deg.list <- NULL
  branch.contrast <- NULL
  ga <- go <- NULL
  if (length(all.branch.ids) == 1) {
    stop(Sys.time(), " [ERROR] There is only one branch in the tree.")
  } else {
    pdata <- object@meta.data[which(object@meta.data$dowsample == 1), c("cell", "branch.id")]
    edata <- object@raw.data[which(object@meta.data$dowsample == 1), ]
    if (is.null(branch.id) & is.null(branch.id.2)) {
      if (verbose) message(Sys.time(), " [INFO] All branches will be calculated.")
      for (bid in all.branch.ids) {
        pdata$contrast <- "go"
        pdata$contrast[which(pdata$branch.id == bid)] = "ga"
        design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
        colnames(design) <- stringr::str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
        fit <- lmFit(t(edata), design)
        contrast <- makeContrasts(ga_go = ga - go,
                                  levels = design)
        fits <- contrasts.fit(fit, contrast)
        ebFit <- eBayes(fits)

        deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = 'fdr', number = Inf)
        deg_sig_list$branch.contrast <- paste0(bid, "_vs_other")
        deg_sig_list$Gene <- rownames(deg_sig_list)
        total.deg.list <- rbind(total.deg.list, deg_sig_list)
      }
    } else if (is.null(branch.id.2)) {
      if (verbose) message(Sys.time(), " [INFO] Some of branches will be calculated.")
      pdata$contrast <- "go"
      pdata$contrast[pdata$branch.id %in% branch.id] = "ga"
      design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
      colnames(design) <- str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
      fit <- lmFit(t(edata), design)
      contrast <- makeContrasts(ga_go = ga - go,
                                levels = design)
      fits <- contrasts.fit(fit, contrast)
      ebFit <- eBayes(fits)

      deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = 'fdr', number = Inf)
      deg_sig_list$branch.contrast <- paste0(paste0(branch.id, collapse = "-"), "_vs_other")
      deg_sig_list$Gene <- rownames(deg_sig_list)
      total.deg.list <- deg_sig_list
    } else {
      if (verbose) message(Sys.time(), " [INFO] Some of branches will be calculated.")
      pdata$contrast <- "gz"
      pdata$contrast[pdata$branch.id %in% branch.id] = "ga"
      pdata$contrast[pdata$branch.id %in% branch.id.2] = "go"
      design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
      colnames(design) <- str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
      fit <- lmFit(t(edata), design)
      contrast <- makeContrasts(ga_go = ga - go,
                                levels = design)
      fits <- contrasts.fit(fit, contrast)
      ebFit <- eBayes(fits)

      deg_sig_list <- topTable(ebFit, coef = 1, adjust.method = 'fdr', number = Inf)
      deg_sig_list$branch.contrast <- paste0(paste0(branch.id, collapse = "-"), "_vs_", paste0(branch.id.2, collapse = "-"))
      deg_sig_list$Gene <- rownames(deg_sig_list)
      total.deg.list <- deg_sig_list
    }
  }
  if (verbose) message(Sys.time(), " [INFO] Calculating differentially expressed markers completed")
  return(total.deg.list)
}
