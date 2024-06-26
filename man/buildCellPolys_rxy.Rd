% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aoo.R
\name{buildCellPolys_rxy}
\alias{buildCellPolys_rxy}
\title{Build simple feature polygons from point data, rotation and shift in x and y direction}
\usage{
buildCellPolys_rxy(thepoints, cellsize, rot, shiftx, shifty)
}
\arguments{
\item{thepoints}{set of points in metres i.e. c(x,y)}

\item{cellsize}{size of cell (length) in metres}

\item{rot}{rotation of the grid in radian}

\item{shiftx}{shift in the x direction in metres}

\item{shifty}{shift in the y direction in metres}
}
\value{
Simple Feature of polygons
}
\description{
Builds cell polygons (as simple features) from points and rotation, shift in X and y returns polygons for ggplot2 and mapping.
Generally used to plot data from AOO calculations.
}
\examples{
#Build and project some points
thepoints <- ptsSquare(19,0.1)
names(thepoints) <- c("lat","long")
thepoints <- simProjWiz(thepoints)
#Check projection information is attributed
attr(thepoints,'crs')
cellsize = 2000
all <- aooFixedRotation(thepoints,cellsize,1296,"ALL")
worstgrid <- all[which.max(all$nofcells),]
worstpoly <- buildCellPolys_rxy(thepoints,cellsize,worstgrid$rotation,worstgrid$xshift,worstgrid$yshift)
ggplot(data=worstpoly) 
  + geom_sf(color="black",fill="red") 
  + geom_point(data=thepoints,aes(X,Y))

}
\references{
Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
}
\author{
Justin Moat. J.Moat@kew.org
}
