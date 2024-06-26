% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ptsSquare}
\alias{ptsSquare}
\title{Builds a set of points in a square area}
\usage{
ptsSquare(noP, gsize = 2000)
}
\arguments{
\item{gsize}{size of square (width or height), (default = 2000)}

\item{nop}{number of points}
}
\value{
dataframe of points (x,y)
}
\description{
Builds a random set of points for a square area
}
\examples{
dfofpoints <- ptsSquare(100,1)
}
\author{
Justin Moat. J.Moat@kew.org
}
