% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sub_loc_frag.R
\name{subLocBuf}
\alias{subLocBuf}
\title{Sub-population or number of locations using the buffer method}
\usage{
subLocBuf(
  thepoints,
  bufferradius = longestAxis(thepoints, returnV = "S")/10,
  returnV = "S"
)
}
\arguments{
\item{thepoints}{dataframe of points of x,y}

\item{bufferradius}{in metres, (default = 1/10 the longest axis of the set of points)}

\item{returnV}{three switches for return values \cr
S = simply the number of sub-pop or locations \cr
AREA = returns the area of the buffered output \cr
SF = return the simple feature for plotting, mapping and export to a GIS}
}
\value{
number of sub-pop/locations
}
\description{
Calculates the number of Sub-population or number of locations using the buffer method
}
\note{
Malin et al 2010 suggest buffer size of 1/10 the longest distance
}
\examples{
#Build some normally distributed point data around the Troodos mountains ~ 10 km diameter
thepoints <- ptsNormal(50,0.1)
#shift to Troodos mountaions
thepoints <- data.frame(long = thepoints$X + 32.8794, lat = thepoints$Y + 34.9220)
#project the points
ppts <- simProjWiz(thepoints)
#get the number of sub-pop/locations using default method
subLocBuf(ppts)
#and area
subLocBuf(ppts,returnV="AREA")
#get the sf object and plot
subpop <- subLocBuf(ppts,returnV="SF")
plot(subpop)
#user defined buffer distance in this case 2 km
subLocBuf(ppts,bufferradius=2000,returnV="S")
}
\references{
Rivers, M. C., Bachman, S. P., Meagher, T. R., Lughadha, E. N., & Brummitt, N. A. (2010). Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data. Biodiversity and Conservation, 19(7), 2071-2085. https://link.springer.com/article/10.1007/s10531-010-9826-9

Willis, F., Moat, J., & Paton, A. (2003). Defining a role for herbarium data in Red List assessments: a case study of Plectranthus from eastern and southern tropical Africa. Biodiversity & Conservation, 12(7), 1537-1552.
}
\seealso{
\code{\link{longestAxis}} to calculate the length of the longest axis
}
\author{
Justin Moat. J.Moat@kew.org
}
