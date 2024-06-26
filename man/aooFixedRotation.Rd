% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aoo.R
\name{aooFixedRotation}
\alias{aooFixedRotation}
\title{Area of Occupancy (AOO) calculated by Systematic shifting and rotating of the grid}
\usage{
aooFixedRotation(
  thepoints,
  cellsize = 2000,
  it = 1296,
  returnV = "S",
  rotation = TRUE
)
}
\arguments{
\item{thepoints}{dataframe of points in metres i.e. c(x,y)}

\item{cellsize}{width of cell in metres (default 2000 m)}

\item{it}{the number of iterations you wish it to run, (default 1296)}

\item{returnV, }{switches to return different sets of results:
S = Simple, returns just the minimum are in km2, (DEFAULT)  
E = Expended simple, returns the solution for the smallest AOO  as list of (area,number of cells, rotation (degrees), shift in x direction, shift in y direction)  
ALL = returns a dataframe of all of the results with (number of cells, rotation (radians), shift in x direction (metres), shift in y direction (metres))  
SF = returns a polygon simple feature for mapping, plotting in ggplot or export to GIS systems.}

\item{rotation}{allow rotation of grids? (default = TRUE). If rotations are selected iterations are shared 50:50 rotation:shift}
}
\value{
dependent on switch, default is area in km2
}
\description{
Calculates the Area of Occupancy AOO (smallest) by shifting and rotating the grid in x and y direction only.
In a very few occasions the minimum solution will not always be achieved but it is quick and consistent (not driven by the number of points).
If your species is near a threshold you may want to increase the number of iterations.
Please cite below if using this algorithm:  
Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.
}
\note{
none
}
\examples{
#Build and project some points
thepoints <- ptsSquare(19,0.1)
names(thepoints) <- c("lat","long")
thepoints <- simProjWiz(thepoints)
attr(thepoints,'crs')
cellsize = 2000

#just get minimum AOO area in km2
aooFixedRotation(thepoints,cellsize,1296,"S")
#get extended results for minimum AOO
aooFixedRotation(thepoints,cellsize,1296,"E")
#get results for all iterations
a <- aooFixedRotation(thepoints,cellsize,1296,"ALL")
hist(a$nofcells)
#Build polygon for plotting
polysrot <- aooFixedRotation(thepoints,cellsize,1296,"SF")
ggplot(data=polysrot) + geom_sf(color="black",fill="red") +
  geom_point(data=thepoints,aes(X,Y))
#Build polygons and plot the worst case
worstgrid <- all[which.max(all$nofcells),]
worstpoly <- buildCellPolys_rxy(thepoints,cellsize,worstgrid$rotation,worstgrid$xshift,worstgrid$yshift)
ggplot(data=worstpoly) + geom_sf(color="black",fill="red") +
  geom_point(data=thepoints,aes(X,Y))

}
\references{
Moat, J., Bachman, S. P., Field, R., & Boyd, D. S. (2018). Refining area of occupancy to address the modifiable areal unit problem in ecology and conservation. Conservation biology, 32(6), 1278-1289.

Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109
}
\seealso{
\code{\link{ratingAoo}} for AOO Ratings

\code{\link{aooFixedGrid}} for fixed grid optimal method

\code{\link{aoo}} for simple AOO method

\code{\link{buildCellPolys_rxy}} for building grid polygons from points, rotation and shift
}
\author{
Justin Moat. J.Moat@kew.org
}
