#' Visualization and analyzation for flow cytometry data
#'
#' Functions and methods to visualize and analyze flow cytometry data.
#'
#' \tabular{ll}{ Package: \tab flowSpy\cr Type: \tab Package\cr Version: \tab
#' 1.2.2\cr Date: \tab 2020-05-07\cr License: \tab GPL-3.0\cr }
#' While high-dimensional single-cell based flow and mass cytometry data has
#' demonstrated increased applications in microenvironment composition and
#' stem-cell research, integrated analyzing workflow design for experimental
#' cytometry data has been challenging. Here, we present flowSpy, an R package
#' designed for the analysis and interpretation of flow and mass cytometry
#' data. We have applied flowSpy to mass cytometry and time course flow
#' cytometry data to validate the usage and practical utility of its
#' computational modules. These use cases introduce flowSpy as a reliable tool
#' for high-dimensional cytometry data workflow and reveal good performance
#' on trajectory reconstruction and pseudotime estimation.
#'
#'
#' @name flowSpy-package
#' @aliases flowSpy-package flowSpy
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
#' vignette(package = "flowSpy")
#' vignette("Quick_start", package = "flowSpy")
#' }
#'
#' @importFrom flowCore read.FCS
#' @importClassesFrom methods ANY character formula logical matrix missing
#' @importMethodsFrom Biobase exprs pData pData<- phenoData sampleNames sampleNames<-
#' @importFrom methods as is
#' @importFrom grid gpar get.gpar grid.points grid.polygon grid.rect current.viewport unit
NULL






