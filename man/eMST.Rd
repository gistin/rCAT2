% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RapoportMP.R
\name{eMST}
\alias{eMST}
\title{Euclidean Minimum spanning tree}
\usage{
eMST(thepoints)
}
\arguments{
\item{thepoints}{dataframe of points in metres i.e. c(X,Y)}
}
\value{
Simple feature of linestring, with a df of X1,Y1,X2,Y2,distance and geom. N.B. X1,Y1 & and X2 Y2 are the to and from points
}
\description{
Calculates the Euclidean Minimum spanning tree from a set of points. 
This is used for the  tree and branch building part of Rapoport's (1982) mean propinquity method.
}
\note{
EMST is computed from Euclidean Minimum Spanning Trees (EMST) using the fast 
Dual-Tree Boruvka algorithm (March, Ram, Gray, 2010, <doi:10.1145/1835804.1835882>) implemented in 'mlpack' - the C++ Machine 
Learning library (Curtin et al., 2013). 'emstreeR' R wrapped for these algorithms  from Allan Quadros & Andre Cancado 2019 
https://cran.r-project.org/web/packages/emstreeR/index.html & https://github.com/allanvc/emstreeR
}
\examples{
#Build and project some points
thepoints <- ptsSquare(19,0.1)
names(thepoints) <- c("lat","long")
thepoints <- simProjWiz(thepoints)
#get the Euclidean Minimum spanning tree
euMST <- eMST(thepoints)
library (ggplot2)
ggplot (data=euMST) + geom_sf(colour="blue")+ geom_point(data=thepoints,aes(X,Y)) 

}
\references{
https://github.com/allanvc/emstreeR

March, W. B., Ram, P., & Gray, A. G. (2010, July). Fast euclidean minimum spanning tree: algorithm, analysis, and applications. In Proceedings of the 16th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 603-612).

Curtin, R. R. et al. (2013). Mlpack: A scalable C++ machine learning library. Journal of Machine Learning Research, v. 14, 2013.

Rivers, M. C., Bachman, S. P., Meagher, T. R., Lughadha, E. N., & Brummitt, N. A. (2010). Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data. Biodiversity and Conservation, 19(7), 2071-2085. https://link.springer.com/article/10.1007/s10531-010-9826-9

Willis, F., Moat, J., & Paton, A. (2003). Defining a role for herbarium data in Red List assessments: a case study of Plectranthus from eastern and southern tropical Africa. Biodiversity & Conservation, 12(7), 1537-1552.

Moat J (2007) Conservation assessment tools extension for ArcView 3.x, version 1.2. GIS Unit, Royal Botanic Gardens, Kew, UK. https://www.kew.org/sites/default/files/2019-02/Conservation_assessment_tools_extension.pdf

Rapoport E.H. 1982. Areography: Geographical Strategies of Species. Pergamon Press, New York.
}
\seealso{
\code{\link{subLocRapoport}} Rapoport's mean propinquity methods
}
\author{
Justin Moat. J.Moat@kew.org
}
