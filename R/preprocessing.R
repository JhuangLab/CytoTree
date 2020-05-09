
###############################################################################
# Functions in preprocessing module are modified from cytofkit package.
# And here is the reference: Hao Chen, Mai Chan Lau, Michael Thomas Wong,
# Evan W. Newell, Michael Poidinger, Jinmiao Chen. Cytofkit: A Bioconductor Package
# for an Integrated Mass Cytometry Data Analysis Pipeline. PLoS Comput Biol, 2016.
###############################################################################

#'
#' Merge the expression matrix from multiple FCS files with preprocessing
#'
#' @name runExprsMerge
#'
#' @description  Apply preprocessing on each FCS file including compensation
#'   (for FCM data only) and transformation with selected markers, then expression
#'   matrix are extracted and merged using one of the methods, \code{all},
#'   \code{min}, \code{fixed} or \code{ceil}
#'
#' @param fcsFiles A vector of FCS file names.
#' @param comp If \verb{TRUE}, does compensation  by compensation matrix contained
#'    in FCS. Agrument also accepts a compensation matrix to be applied.
#'    Otherwise \verb{FALSE}.
#' @param transformMethod Data Transformation method, including \code{autoLgcl},
#'    \code{cytofAsinh}, \code{logicle} and \code{arcsinh}, or \code{none} to
#'    avoid transformation.
#' @param scaleTo Scale the expression to a specified range c(a, b), default is NULL.
#' @param mergeMethod Merge method for mutiple FCS expression data. cells can be
#'    combined using one of the four different methods including \code{ceil},
#'    \code{all}, \code{min}, \code{fixed}. The default option is \code{ceil},
#'    up to a fixed number (specified by \code{fixedNum}) of cells are sampled
#'    without replacement from each fcs file and combined for analysis.
#'    \code{all}: all cells from each fcs file are combined for analysis.
#'    \code{min}: The minimum number of cells among all the selected fcs files
#'    are sampled from each fcs file and combined for analysis.
#'    \code{fixed}: a fixed num (specified by fixedNum) of cells are sampled
#'    (with replacement when the total number of cell is less than fixedNum)
#'    from each fcs file and combined for analysis.
#' @param fixedNum The fixed number of cells to be extracted from each FCS file.
#' @param ... Other arguments passed to \code{runExprsExtract}
#'
#' @return A matrix containing the merged expression data, with selected markers.
#' @seealso \code{\link{runExprsExtract}}
#' @export
#'
#' @author Chen Hao
#' @references Hao Chen, Mai Chan Lau, Michael Thomas Wong, Evan W. Newell,
#'    Michael Poidinger, Jinmiao Chen. Cytofkit: A Bioconductor Package for
#'    an Integrated Mass Cytometry Data Analysis Pipeline. PLoS Comput Biol, 2016.
#'
#' @examples
#'
#'
#' if (FALSE) {
#'   # See vignette tutorials for more information
#'   vignette("Quick_start", package = "flowSpy")
#'
#'   # Path to your FCS files
#'   fcs.path <- "flowSpy-dataset/FCS/usecase2/"
#'   fcs.files <- paste0(fcs.path, "D", c(0,2,4,6,8,10), "-sub.fcs")
#'
#'   # Merge FCS files, and each file contain 2000 cells
#'   set.seed(1)
#'   fcs.data <- runExprsMerge(fcs.files, comp = F, transformMethod = "none", fixedNum = 2000)
#' }
#'
#'
runExprsMerge <- function(fcsFiles,
                          comp = FALSE,
                          transformMethod = c("autoLgcl", "cytofAsinh", "logicle", "arcsinh", "logAbs", "none"),
                          scaleTo = NULL,
                          mergeMethod = c("ceil", "all", "fixed", "min"),
                          fixedNum = 2000, ...) {

  transformMethod <- match.arg(transformMethod)
  mergeMethod <- match.arg(mergeMethod)

  exprsL <- mapply(runExprsExtract, fcsFiles,
                   MoreArgs = list(comp = comp,
                                   transformMethod = transformMethod,
                                   scaleTo = scaleTo, ...),
                   SIMPLIFY = FALSE)

  ## test if number of events in any fcs less than fixedNum
  eventCountTest <- suppressWarnings(any(lapply(exprsL, function(x) if (nrow(x) < fixedNum) {1} else {0})))
  ## solution 1: change mergeMethod from fixed to ceil
  #if(mergeMethod == "fixed" && eventCountTest == TRUE){
  #  mergeMethod <- "ceil"
  #}
  ## solution 2: use lowest number of fcs events as fixedNum
  if(mergeMethod == "fixed" && eventCountTest == TRUE){
    warning("One or more FCS files have less events than specified fixedNum; using lowest as fixedNum")
    fixedNum <- min(rapply(exprsL, nrow))
  }
  switch(mergeMethod,
         ceil = {
           mergeFunc <- function(x) {
             if (nrow(x) < fixedNum) {
               x
             } else {
               x[sample(nrow(x), size = fixedNum, replace = FALSE), , drop = FALSE]
             }
           }
           merged <- do.call(rbind, lapply(exprsL, mergeFunc))
         },
         all = {
           merged <- do.call(rbind, exprsL)
         },
         fixed = {
           mergeFunc <- function(x) {
             x[sample(nrow(x), size = fixedNum, replace = ifelse(nrow(x) < fixedNum, TRUE, FALSE)), , drop=FALSE]
           }
           merged <- do.call(rbind, lapply(exprsL, mergeFunc))
         },
         min = {
           minSize <- min(sapply(exprsL, nrow))
           mergeFunc <- function(x) {
             x[sample(nrow(x), size = minSize, replace = FALSE), , drop=FALSE]
           }
           merged <- do.call(rbind, lapply(exprsL, mergeFunc))
         })

  return(merged)
}




#'
#' Extract the expression data from a FCS file with preprocessing
#'
#' @name runExprsExtract
#'
#' @description  Extract the FCS expresssion data with preprocessing of
#'    compensation (for FCM data only) and transformation. Transformtion
#'    methods includes \code{autoLgcl}, \code{cytofAsinh}, \code{logicle}
#'    (customizable) and \code{arcsinh} (customizable).
#'
#' @param fcsFile The name of the FCS file.
#' @param verbose If \verb{TRUE}, print the message details of FCS loading.
#' @param comp If \verb{TRUE}, does compensation  by compensation matrix
#'    contained in FCS. Agrument also accepts a compensation matrix to be
#'    applied. Otherwise \verb{FALSE}.
#' @param transformMethod Data Transformation method, including \code{autoLgcl},
#'    \code{cytofAsinh}, \code{logicle} and \code{arcsinh}, or \code{none}
#'    to avoid transformation.
#' @param showDesc logical. Whether to show \code{desc} name in the output matrix.
#' @param keepRaw logical. Whether to keep raw data for FSC and SSC.
#' @param scaleTo Scale the expression to a specified range c(a, b), default is NULL.
#' @param q Quantile of negative values removed for auto w estimation,
#'    default is 0.05, parameter for autoLgcl transformation.
#' @param l_w Linearization width in asymptotic decades, parameter for
#'    logicle transformation.
#' @param l_t Top of the scale data value, parameter for logicle transformation.
#' @param l_m Full width of the transformed display in asymptotic decades,
#'    parameter for logicle transformation.
#' @param l_a Additional negative range to be included in the display
#'    in asymptotic decades, parameter for logicle transformation.
#' @param a_a Positive double that corresponds to the base of the arcsinh
#'    transformation, \code{arcsinh} = asinh(a + b * x) + c).
#' @param a_b Positive double that corresponds to a scale factor of the
#'    arcsinh transformation, \code{arcsinh} = asinh(a + b * x) + c).
#' @param a_c Positive double that corresponds to another scale factor
#'    of the arcsinh transformation, \code{arcsinh} = asinh(a + b * x) + c).
#'
#' @return A transformed expression data matrix
#'
#' @importFrom flowCore read.FCS compensate estimateLogicle logicleTransform
#' @importFrom flowCore parameters transformList arcsinhTransform biexponentialTransform
#' @importMethodsFrom flowCore transform
#' @importClassesFrom flowCore transformList
#'
#' @export
#'
#' @author Chen Hao
#' @references Hao Chen, Mai Chan Lau, Michael Thomas Wong, Evan W. Newell,
#'    Michael Poidinger, Jinmiao Chen. Cytofkit: A Bioconductor Package for
#'    an Integrated Mass Cytometry Data Analysis Pipeline. PLoS Comput Biol, 2016.
#'
#' @examples
#'
#' if (FALSE) {
#' # See vignette tutorials for more information
#' vignette(package = "flowSpy")
#' vignette("Quick_start", package = "flowSpy")
#'
#'   # Path to your FCS files
#'   fcs.path <- "flowSpy-dataset/FCS/usecase1/"
#'   fcs.file <- paste0(fcs.path, "FR-FCM-ZY9R-Bone_Marrow_cytof.fcs")
#'
#'   # Read FCS files
#'   exp.data <- runExprsExtract(fcs.file, showDesc = FALSE, transformMethod = "autoLgcl")
#' }
#'
#'
runExprsExtract <- function(fcsFile,
                            verbose = FALSE,
                            comp = FALSE,
                            transformMethod = c("autoLgcl", "cytofAsinh", "logicle", "arcsinh", "logAbs", "none"),
                            scaleTo = NULL,
                            showDesc = TRUE,
                            keepRaw = TRUE,
                            q = 0.05,
                            l_w = 0.1, l_t = 4000, l_m = 4.5, l_a = 0,
                            a_a = 1, a_b = 1, a_c =0) {

  transformMethod <- match.arg(transformMethod)

  ## load FCS files
  name <- sub(".fcs$", "", basename(fcsFile))
  if (verbose) {
    fcs <- read.FCS(fcsFile, transformation = FALSE)
  } else {
    fcs <- suppressWarnings(read.FCS(fcsFile, transformation = FALSE))
  }

  ## compensation
  if(is.matrix(comp) || is.data.frame(comp)){
    fcs <- applyComp(fcs, comp)
    message("    Compensation is applied on", fcsFile, "\n")
  }else if(isTRUE(comp)) {
    if(!is.null(fcs@description$SPILL)) {
      fcs <- applyComp(fcs, fcs@description[["SPILL"]])
      message("    Compensation is applied on ", fcsFile, "\n")
    }else if(!is.null(fcs@description$SPILLOVER)) {
      fcs <- applyComp(fcs, fcs@description[["SPILLOVER"]])
      message("    Compensation is applied on ", fcsFile, "\n")
    }else if(!is.null(fcs@description$COMP)) {
      fcs <- applyComp(fcs, fcs@description[["COMP"]])
      message("    Compensation is applied on ", fcsFile, "\n")
    }else{
      warning(Sys.time(), " [WARNING] ", "Cannot find compensation matrix in the FCS files!
              Please CHECK the keyword of 'SPILL', 'SPILLOVER', or 'COMP'
              in the FCS file and make sure it stores the compensation matrix.")
    }
  }

  ## match marker names to get marker ID, use all if NULL
  pd <- fcs@parameters@data

  ## Exclude "Time", "Event" channel
  exclude_channels <- grep("Time|Event", colnames(fcs@exprs), ignore.case = TRUE)
  marker_id <- setdiff(seq_along(colnames(fcs@exprs)), exclude_channels)

  size_channels <- grep("FSC|SSC", colnames(fcs@exprs), ignore.case = TRUE)
  transMarker_id <- setdiff(marker_id, size_channels)

  ## exprs transformation
  switch(transformMethod,
         cytofAsinh = {
           data <- fcs@exprs
           data[ ,transMarker_id] <- apply(data[ ,transMarker_id, drop=FALSE], 2, cytofAsinh)
           exprs <- data[ ,marker_id, drop=FALSE]
         },
         autoLgcl = {
           trans <- autoLgcl(fcs, channels = colnames(fcs@exprs)[transMarker_id], q = q)
           transformed <- flowCore::transform(fcs, trans)
           exprs <- transformed@exprs[, marker_id, drop=FALSE]
         },
         logicle = {
           data <- fcs@exprs
           trans <- flowCore::logicleTransform(w = l_w, t = l_t, m = l_m, a = l_a)
           data[ ,transMarker_id] <- apply(data[ ,transMarker_id, drop=FALSE], 2, trans)
           exprs <- data[ ,marker_id, drop=FALSE]
         },
         arcsinh = {
           data <- fcs@exprs
           trans <- flowCore::arcsinhTransform(a = a_a, b = a_b, c = a_c)
           data[ ,transMarker_id] <- apply(data[ ,transMarker_id, drop=FALSE], 2, trans)
           exprs <- data[ ,marker_id, drop=FALSE]
         },
         logAbs = {
           data <- fcs@exprs
           cs <- apply(abs(data[ ,transMarker_id]), 2, sum)
           transMarker_id <- transMarker_id[cs > 0]
           cs <- cs[cs > 0]
           norm_factors <- (10**ceiling(log10(median(cs))))/cs
           data[ ,transMarker_id] <- round(log10(sweep(abs(data[ ,transMarker_id]), 2, norm_factors, "*")+1), digits=3)
           exprs <- data[ ,marker_id, drop=FALSE]
         },
         none = {
           data <- fcs@exprs
           exprs <- data[ ,marker_id, drop=FALSE]
         })

  ## apply linear transformation for the "FSC-x", "SSC-x" channel if exists
  if((length(size_channels) > 0) & (keepRaw == FALSE)){
    if(any(size_channels %in% marker_id)){
      used_size_channel <- size_channels[size_channels %in% marker_id]
      used_size_channel_id <- match(used_size_channel, marker_id)
      exprs[ ,used_size_channel_id] <- apply(exprs[ , used_size_channel_id, drop=FALSE], 2,
                                             function(x) scaleData(x, range=c(0, 4.5)))
    }
  }

  ## rescale all data to same range
  if (!is.null(scaleTo)) {
    exprs <- apply(exprs, 2, function(x) scaleData(x, scaleTo))
  }

  ## add rownames and colnames
  if (showDesc) {
    col_names <- paste0(pd$name, "<", pd$desc,">")
  } else {
    col_names <- paste0(pd$name)
  }

  colnames(exprs) <- col_names[marker_id]
  row.names(exprs) <- paste(name, 1:nrow(exprs), sep = "_")

  return(exprs)
}


#' apply compensation on the FCS expression data
#'
#' @param fcs FCS file.
#' @param compMatrix Compensation matrix.
#' @noRd
#' @return Compensated expression value
applyComp <- function(fcs, compMatrix) {
  comp_fcs <- compensate(fcs, compMatrix)
}

#' rescale the data
#'
#' @param x data.
#' @param range The range of the data.
#' @noRd
#' @return scaled data
#'
#' @author Chen Hao
#' @references Hao Chen, Mai Chan Lau, Michael Thomas Wong, Evan W. Newell,
#'    Michael Poidinger, Jinmiao Chen. Cytofkit: A Bioconductor Package for
#'    an Integrated Mass Cytometry Data Analysis Pipeline. PLoS Comput Biol, 2016.
#'
scaleData <- function(x, range = c(0, 4.5)) {
  (x - min(x))/(max(x) - min(x)) * (range[2] - range[1]) + range[1]
}

#'
#' Noise reduced arsinh transformation
#'
#' Inverse hyperbolic sine transformation (arsinh) with a cofactor of 5,
#' reduce noise from negative values. Adopted from Plos Comp reviewer
#'
#' @param value A vector of numeric values.
#' @param cofactor Cofactor for asinh transformation, default 5 for mass cytometry data.
#' @noRd
#' @return transformed value
#' @importFrom stats rnorm median quantile
#'
#' @author Chen Hao
#' @references Hao Chen, Mai Chan Lau, Michael Thomas Wong, Evan W. Newell,
#'    Michael Poidinger, Jinmiao Chen. Cytofkit: A Bioconductor Package for
#'    an Integrated Mass Cytometry Data Analysis Pipeline. PLoS Comput Biol, 2016.
#'
#'
cytofAsinh <- function(value, cofactor = 5) {
  value <- value-1
  loID <- which(value < 0)
  if(length(loID) > 0)
    value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
  value <- value / cofactor
  value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
  return(value)
}


#' a modified version of "estimateLogicle" from flowCore
#'
#' Used boxplot outlier detection to filter outliers in negative values
#' before calculating the r using the fifth percentile of the negative values.
#'
#' @param x A flowFrame object.
#' @param channels Channel names to be transformed.
#' @param m The full width of the transformed display in asymptotic decades.
#'    \code{m} should be greater than zero.
#' @param q The percentile of negative values used as reference poiont of negative range.
#' @importFrom methods is
#' @importFrom flowCore logicleTransform
#' @importFrom stats IQR
#' @noRd
#' @return a list of autoLgcl transformations
#'
#' @author Chen Hao
#' @references Hao Chen, Mai Chan Lau, Michael Thomas Wong, Evan W. Newell,
#'    Michael Poidinger, Jinmiao Chen. Cytofkit: A Bioconductor Package for
#'    an Integrated Mass Cytometry Data Analysis Pipeline. PLoS Comput Biol, 2016.
#'
autoLgcl <- function(x, channels, m = 4.5, q = 0.05) {
  if (!is(x, "flowFrame"))
    stop("x has to be an object of class \"flowFrame\"")
  if (missing(channels))
    stop("Please specify the channels to be logicle transformed")
  indx <- channels %in% colnames(x@exprs)
  if (!all(indx))
    stop(paste("Channels", channels[!indx], "were not found in the FCS file.\n ",
               sep = " "))

  trans <- lapply(channels, function(p) {
    data <- x@exprs[, p]
    w <- 0
    t <- max(data)
    ndata <- data[data < 0]
    ## use 1.5 * IQR to filter outliers in negative values
    nThres <- quantile(ndata, 0.25) - 1.5 * stats::IQR(ndata)
    ndata <- ndata[ndata >= nThres]
    transId <- paste(p, "autolgclTransform", sep = "_")

    if (length(ndata)) {
      r <- .Machine$double.eps + quantile(ndata, q)
      ## Check to avoid failure of negative w
      if (10^m * abs(r) <= t) {
        w <- 0
      } else {
        w <- (m - log10(t/abs(r)))/2
        if(is.nan(w) || w>2) {
          warning(Sys.time(), " [WARNING] ", paste0("autoLgcl failed for channel: ", p, "; using default logicle transformation!"))
          w <- 0.1
          t <- 4000
          m <- 4.5
        }
      }
    }
    logicleTransform(transformationId = transId,
                     w = w, t = t, m = m, a = 0)
  })
  transformList(channels, trans)
}



