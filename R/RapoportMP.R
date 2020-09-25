#Scripts for Rapoport's mean propinquity

###########################################################
#returns the Euclidean Minimum spanning tree from a set of points
###########################################################
#' @title Euclidean Minimum spanning tree 
#' @description
#' Calculates the Euclidean Minimum spanning tree from a set of points. 
#' This is used for the  tree and branch building part of Rapoport's (1982) mean propinquity method. 
#' 
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points in metres i.e. c(X,Y)
#' @return Simple feature of linestring, with a df of X1,Y1,X2,Y2,distance and geom. N.B. X1,Y1 & and X2 Y2 are the to and from points 
#' @note 
#' EMST is computed from Euclidean Minimum Spanning Trees (EMST) using the fast 
#' Dual-Tree Boruvka algorithm (March, Ram, Gray, 2010, <doi:10.1145/1835804.1835882>) implemented in 'mlpack' - the C++ Machine 
#' Learning library (Curtin et al., 2013). 'emstreeR' R wrapped for these algorithms  from Allan Quadros & Andre Cancado 2019 
#' https://cran.r-project.org/web/packages/emstreeR/index.html & https://github.com/allanvc/emstreeR
#' @examples
#'#Build and project some points
#'thepoints <- ptsSquare(19,0.1)
#'names(thepoints) <- c("lat","long")
#'thepoints <- simProjWiz(thepoints)
#'#get the Euclidean Minimum spanning tree
#'euMST <- eMST(thepoints)
#'library (ggplot2)
#'ggplot (data=euMST) + geom_sf(colour="blue")+ geom_point(data=thepoints,aes(X,Y)) 
#'
#' @export
#' @import sf
#' @import emstreeR
#' @references
#' https://github.com/allanvc/emstreeR
#' 
#' March, W. B., Ram, P., & Gray, A. G. (2010, July). Fast euclidean minimum spanning tree: algorithm, analysis, and applications. In Proceedings of the 16th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 603-612).
#' 
#' Curtin, R. R. et al. (2013). Mlpack: A scalable C++ machine learning library. Journal of Machine Learning Research, v. 14, 2013.
#' 
#' Rivers, M. C., Bachman, S. P., Meagher, T. R., Lughadha, E. N., & Brummitt, N. A. (2010). Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data. Biodiversity and Conservation, 19(7), 2071-2085. https://link.springer.com/article/10.1007/s10531-010-9826-9
#' 
#' Willis, F., Moat, J., & Paton, A. (2003). Defining a role for herbarium data in Red List assessments: a case study of Plectranthus from eastern and southern tropical Africa. Biodiversity & Conservation, 12(7), 1537-1552.
#' 
#' Moat J (2007) Conservation assessment tools extension for ArcView 3.x, version 1.2. GIS Unit, Royal Botanic Gardens, Kew, UK. https://www.kew.org/sites/default/files/2019-02/Conservation_assessment_tools_extension.pdf
#' 
#' Rapoport E.H. 1982. Areography: Geographical Strategies of Species. Pergamon Press, New York.
#' @seealso \code{\link{subLocRapoport}} Rapoport's mean propinquity methods

eMST <- function (thepoints){
  edges <- nrow(thepoints)- 1
  cmst <- ComputeMST (thepoints,verbose = FALSE)[1:edges,3:5] #drops the last point and the X/Y's
  #build the simple feature, probably better ways of doing this, but it works
  from <- thepoints[cmst$from,]
  to <- thepoints[cmst$to,]
  df <- cbind(from,to,cmst$distance)
  names(df) <- c("X1","Y1","X2","Y2","distance")
  df$geom <- st_sfc(sapply(1:nrow(df),function(i){st_linestring(t(matrix(unlist(df[i,1:4]), 2, 2)))},simplify = FALSE))
  sfmst <- st_sf(df)
  st_crs(sfmst) <- attr(thepoints,'crs')
  sfmst
}


###########################################################
#returns Rapoport's mean propinquity area or number of sub-populations from a set of points
###########################################################
#' @title Sub-population or number of locations using Rapoport's mean propinquity
#' @description
#' Calculates Rapoport's (1982) mean propinquity from a set of points. 
#' For details on Rapoport’s mean propinquity methodsee Willis et al. 2003 and Rapoport 1982. This technique is
#' based on a Euclidean minimum spanning tree (eMST), which is a
#' set of lines (branches) that connects all points in the minimum possible distance. Sub-populations
#' are defined when the branches (line from a pair of points) distance is greater than the defined distance 
#' (Rapoport suggests twice the mean edge distance (Willis et al. 2003)). 
#' Area is defined by  buffering this tree by a defined distance (Rapoport suggests the mean) 
#' any isolated populations (single populations or points) are also buffered by this distance.
#' 
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points in metres i.e. c(X,Y)
#' @param barrierDis distance (in m) to be used to define barriers (if not defined the default is 2 x mean branch length of the Euclidean Minimum spanning tree). Rivers et al (2010) and Willis et al (2003) suggest 1/10 of the longest Axis may be useful.
#' @param bufferDis distance (in m) to be used to buffer points and connected edges (if not defined the default is the mean branch length of the Euclidean Minimum spanning tree )
#' @param returnV, switches to return different sets of results: \cr
#' S = Simple, returns just the number of sub-populations or locations (default) \cr
#' AREA = returns the area (km2) of the buffered points and tree \cr
#' SF = returns a list with two simple features (tree and buffers) for visualisation, mapping, further analysis and export to a GIS
#' @return dependent on switch, default is area in km2
#' @return see returnV
#' @examples
#'#Build and project some points
#'thepoints <- ptsSquare(19,0.1)
#'names(thepoints) <- c("lat","long")
#'thepoints <- simProjWiz(thepoints)
#'#number of sub-populations/locations using default distances
#'subLocRapoport(thepoints)
#'#defining your own distances
#'subLocRapoport(thepoints,barrierDis = 2000,bufferDis = 1000)
#'#area of sub populations/locations using default distances
#'subLocRapoport(thepoints,returnV = "AREA")
#'#results as simple features for plotting
#'sfs <- subLocRapoport(thepoints,returnV = "SF")
#'library(ggplot2)
#'ggplot (data=sfs$tree) + geom_sf(aes(color=barrier)) + geom_sf(data=sfs$buffers,colour='green', fill=NA) + geom_point(data=thepoints,aes(X,Y)) 
#'#get the mean of distance for the tree branches
#'mean(sfs$tree$distance)
#'
#' @export
#' @import sf
#' @import emstreeR
#' @references
#' https://github.com/allanvc/emstreeR
#' 
#' March, W. B., Ram, P., & Gray, A. G. (2010, July). Fast euclidean minimum spanning tree: algorithm, analysis, and applications. In Proceedings of the 16th ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 603-612).
#' 
#' Rivers, M. C., Bachman, S. P., Meagher, T. R., Lughadha, E. N., & Brummitt, N. A. (2010). Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data. Biodiversity and Conservation, 19(7), 2071-2085. https://link.springer.com/article/10.1007/s10531-010-9826-9
#' 
#' Willis, F., Moat, J., & Paton, A. (2003). Defining a role for herbarium data in Red List assessments: a case study of Plectranthus from eastern and southern tropical Africa. Biodiversity & Conservation, 12(7), 1537-1552.
#' 
#' Moat J (2007) Conservation assessment tools extension for ArcView 3.x, version 1.2. GIS Unit, Royal Botanic Gardens, Kew, UK. https://www.kew.org/sites/default/files/2019-02/Conservation_assessment_tools_extension.pdf
#' 
#' Rapoport E.H. 1982. Areography: Geographical Strategies of Species. Pergamon Press, New York. 
#' @seealso \code{\link{eMST}} Euclidean Minimum spanning tree 
#' @seealso \code{\link{longestAxis}} Longest distance from a set of points 


subLocRapoport <- function(thepoints,barrierDis,bufferDis,returnV='S'){
  #get the eMST
  euMST <- eMST(thepoints)
  #check and build defaults for distance measures
  if (missing(barrierDis)){barrierDis <- mean(euMST$distance)*2}
  if (missing(bufferDis)){bufferDis <- mean(euMST$distance)}
  #mark the barriers
  euMST$barrier <- euMST$distance > barrierDis
  #buffer branches
  branch_area <- st_union(st_buffer(euMST[!euMST$barrier,],bufferDis))
  #buffer points
  ppts_sf <- st_as_sf(thepoints, coords = c("X", "Y"), crs = attr(thepoints,'crs'))
  point_area <- st_union(st_buffer(ppts_sf,bufferDis))
  pop_poly <- st_union(branch_area,point_area)
  if (returnV=='SF'){return(list(tree=euMST,buffers=pop_poly))}
  if (returnV=='AREA'){
    parea <- st_area(pop_poly)
    units(parea) <- "km^2"
    return(parea)
    }
  return(length(st_geometry(pop_poly)[[1]]))
}