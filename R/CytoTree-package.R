#' Visualization and analyzation for flow cytometry data
#'
#' Functions and methods to visualize and analyze flow cytometry data.
#'
#' \tabular{ll}{ Package: \tab CytoTree \cr Type: \tab Package \cr Version: \tab
#' 1.1.2 \cr Date: \tab 2020-11-01 \cr License: \tab GPL-3.0 \cr }
#' While high-dimensional single-cell based flow and mass cytometry data has
#' demonstrated increased applications in microenvironment composition and
#' stem-cell research, integrated analyzing workflow design for experimental
#' cytometry data has been challenging. Here, we present CytoTree, an R package
#' designed for the analysis and interpretation of flow and mass cytometry
#' data. We have applied CytoTree to mass cytometry and time course flow
#' cytometry data to validate the usage and practical utility of its
#' computational modules. These use cases introduce CytoTree as a reliable tool
#' for high-dimensional cytometry data workflow and reveal good performance
#' on trajectory reconstruction and pseudotime estimation.
#'
#'
#' @name CytoTree-package
#' @aliases CytoTree-package CytoTree
#' @docType package
#' @author
#' Maintainer: Yuting Dai <forlynna@@sjtu.edu.cn>
#' Authors: Yuting Dai
#' @keywords package
#' @examples
#'
#' if (FALSE) {
#' ## examples go here
#' ## See vignette tutorials
#' vignette(package = "CytoTree")
#' }
#'
#' @importFrom flowCore read.FCS
#' @importClassesFrom methods ANY character formula logical matrix missing
#' @importMethodsFrom Biobase exprs pData pData<- phenoData sampleNames sampleNames<-
#' @importFrom methods as is
#' @importFrom grid gpar get.gpar grid.points grid.polygon grid.rect current.viewport unit
NULL






