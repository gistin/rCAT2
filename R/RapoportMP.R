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
#'
eMST <- function (thepoints){
  crs <- attr(thepoints, "crs")
  
  edges <- nrow(thepoints) - 1
  cmst <- ComputeMST(thepoints, verbose=FALSE)
  cmst <- cmst[1:edges, 3:5] #drops the last point and the X/Y's
  
  # not sure this is the best way, but hopefully clear
  line_idx <- split(cmst, f=1:nrow(cmst))
  
  lines <- lapply(line_idx, function(idx) {
    p1 <- c(thepoints[idx$from,]$X, thepoints[idx$from,]$Y)
    p2 <- c(thepoints[idx$to,]$X, thepoints[idx$to,]$Y)
    
    st_linestring(rbind(p1, p2))
  })
  
  lines <- st_sfc(lines, crs=crs)
  
  st_sf(
    X1=thepoints[cmst$from,]$X,
    Y1=thepoints[cmst$from,]$Y,
    X2=thepoints[cmst$to,]$X,
    Y2=thepoints[cmst$to,]$Y,
    distance=cmst$distance,
    geometry=lines
  )
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
#'
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
  
  if (returnV == 'SF') {
    list(tree=euMST, buffers=pop_poly)
  } else if (returnV == 'AREA') {
    parea <- st_area(pop_poly)
    
    if (is.null(attr(area, "units"))) {
      area <- area / 1e6
    }
    
    units(parea) <- "km^2"
    parea
  } else {
    polygons <- st_cast(pop_poly, "POLYGON")
    length(polygons)
  }
}
