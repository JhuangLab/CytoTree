% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{fetchCell}
\alias{fetchCell}
\title{Fetching cellls of CYT}
\usage{
fetchCell(object, logical.connect = "or", verbose = FALSE, ...)
}
\arguments{
\item{object}{An CYT object}

\item{logical.connect}{character. "and" or "or"}

\item{verbose}{logical. Whether to print calculation progress.}

\item{...}{Paramters to pass to limitation}
}
\value{
a vector containing cell names
}
\description{
Fetching cellls of CYT
}
\examples{

cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
cyt <- readRDS(file = cyt.file)

cell.fetch <- fetchCell(cyt, traj.value.log = 0.01)
cell.fetch <- fetchCell(cyt, stage = c("D0", "D10"))
cell.fetch <- fetchCell(cyt, stage = c("D0", "D10"), traj.value.log = 0.01,
                        logical.connect = "or")


}
