% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.R
\name{runMclust}
\alias{runMclust}
\title{runMclust}
\usage{
runMclust(object, scale = FALSE, verbose = FALSE, ...)
}
\arguments{
\item{object}{a CYT object}

\item{scale}{logical. Whether to use scaled data in Mclust.}

\item{verbose}{logical. Whether to print calculation progress.}

\item{...}{Parameters passing to \code{\link[mclust]{Mclust}} function}
}
\value{
a CYT object with mclust.id in meta.data
}
\description{
Model-based clustering based on parameterized finite Gaussian mixture models.
   This function is based on \code{\link[mclust]{Mclust}}.
}
\examples{

cyt.file <- system.file("extdata/cyt.rds", package = "CytoTree")
cyt <- readRDS(file = cyt.file)

cyt <- runMclust(cyt, verbose = TRUE)


}
\seealso{
\code{\link[mclust]{Mclust}}
}
