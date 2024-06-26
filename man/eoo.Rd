% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eoo.R
\name{eoo}
\alias{eoo}
\title{Extent of Occurrence (EOO) Area}
\usage{
eoo(points, returnV = "S")
}
\arguments{
\item{returnV}{switch to return different sets of results: \cr
S = Simple, returns just the minimum area in km2, (DEFAULT) \cr
SF = returns a polygon simple feature for mapping, plotting in ggplot/plot or to export to a GIS format}

\item{thepoints}{dataframe of points in metres i.e. c(x,y)}
}
\value{
float_value area of EOO polygon or sf polygon
}
\description{
Calculates the Extent of Occurrence in km2 or returns a simple feature polygon from a set of points (x,y)
}
\note{
area returned is in x,y units, but negative as polygon is constructed anticlockwise
}
\examples{
x <- runif (20,0,10)
y <- runif (20,0,10)
df <- data.frame(x,y) 
eoo (df)
#######
spoly <- eoo (df,TRUE)
plot(spoly)
}
\references{
Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 

Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591
}
\seealso{
\code{\link{ratingEoo}} for EOO Ratings
}
\author{
Justin Moat. J.Moat@kew.org
}
