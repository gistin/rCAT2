% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aoo.R
\name{aoo}
\alias{aoo}
\title{Area of Occupancy (AOO), grid orgin 0,0}
\usage{
aoo(thepoints, cellsize = 2000, returnV = "S")
}
\arguments{
\item{thepoints}{set of points in metres i.e. c(x,y)}

\item{cellsize}{size of cell (length) in metres}

\item{returnV, }{switches to return different sets of results: \cr
S = Simple, returns just the AOO area in km2, (DEFAULT) \cr
E = Expended simple, returns the solution for the AOO as list of (area,number of cells, rotation (0 degrees), shift in x direction(0), shift in y direction(0)). This is used so as be compatiable with other AOO calculations. \cr
SF = returns a polygon simple feature for mapping and plotting in ggplot/plot or export to GIS format.}
}
\value{
as returnV, default is area in km2
}
\description{
Calculates the number area the of occupied cells for (Area of Occupancy AOO) from a set of points (x,y), projected into metres, with origin 0,0. 
Please cite: Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289. if using this algorithm:
}
\details{
calculates a very simple AOO area from a set of points
}
\examples{
library(ggplot2)
#Build and project some points
thepoints <- ptsSquare(19,0.1)
names(thepoints) <- c("lat","long")
thepoints <- simProjWiz(thepoints)
attr(thepoints,'crs')
cellsize = 2000

#return area in km2
aoo (thepoints,cellsize)
#return list of parameters
aoo (thepoints,cellsize,returnV="E")
#return polygon for plotting
gridpoly <- aoo(thepoints,cellsize,returnV="SF")
ggplot(data=gridpoly) 
  + geom_sf(color="black",fill="red") 
  + geom_point(data=thepoints,aes(X,Y))
}
\references{
Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.

Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., (2011). Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109
}
\seealso{
\code{\link{ratingAoo}} for AOO Ratings from IUCN categories
}
\author{
Justin Moat. J.Moat@kew.org
}
