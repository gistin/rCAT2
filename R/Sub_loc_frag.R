#Scripts for Sub-population and number of locations (same methods) or fragmentation
##################################################################################
#calculates the number of Sub-population or number of locations from buffer method
##################################################################################
#' @title Sub-population or number of locations using buffer method
#' @description 
#' Calculates the number of Sub-population or number of locations from buffer method
#' @author Justin Moat. J.Moat@kew.org
#' @note Malin et al 2010 suggest buffer size of 1/10 the longest distance
#' @param thepoints dataframe of points of x,y
#' @param bufferradius in metres, (default = 1/10 the longest axis of the set of points)
#' @param returnV three switches for return values 
#' S = simply the number of sub-pop or locations 
#' AREA = returns the area of the buffered output
#' SF = return the simple feature for plotting and mapping
#' 
#' @return number of sub-pop/locations
#' @examples 
#'#Build some normally distributed point data around the Troodos mountains ~ 10 km diameter
#'thepoints <- normalofPs(50,0.1)
#'#shift to Troodos mountaions
#'thepoints <- data.frame(long = thepoints$X + 32.8794, lat = thepoints$Y + 34.9220)
#'#project the points
#'ppts <- simProjWiz(thepoints)
#'#get the number of sub-pop/locations using default method
#'subLocBuf(ppts)
#'#and area
#'subLocBuf(ppts,returnV="AREA")
#'#get the sf object and plot
#'subpop <- subLocBuf(ppts,returnV="SF")
#'plot(subpop)
#'#user defined buffer distance in this case 2 km
#'subLocBuf(ppts,bufferradius=2000,returnV="S")
#' @export
#' @import sf
#' @references 
#' Rivers, M. C., Bachman, S. P., Meagher, T. R., Lughadha, E. N., & Brummitt, N. A. (2010). Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data. Biodiversity and Conservation, 19(7), 2071-2085. https://link.springer.com/article/10.1007/s10531-010-9826-9
#' 
#' Willis, F., Moat, J., & Paton, A. (2003). Defining a role for herbarium data in Red List assessments: a case study of Plectranthus from eastern and southern tropical Africa. Biodiversity & Conservation, 12(7), 1537-1552.
#' @seealso \code{\link{longestAxis}} to calculate the length of the longest axis

subLocBuf <- function (thepoints,bufferradius=longestAxis(thepoints,returnV='S')/10,returnV="S"){
  thepointsSF <- st_as_sf(thepoints, coords = c("X", "Y"),crs= attr(thepoints,'crs'))
  buf <- st_union(st_buffer(thepointsSF,bufferradius))
  polys <- st_cast(buf,"POLYGON")
  if(returnV=="SF"){return (polys)}
  if(returnV=="AREA"){
    parea <- sum(st_area(polys))
    units(parea) <- "km^2"
    return(parea)
  }
  #message(bufferradius)
  length(polys)
}

##################################################################################
#calculates the number of Sub-population or number of locations from cell adjacency method
##################################################################################
#' @title Sub-population or number of locations using grid adjacency
#' @description 
#' Calculates the number of Sub-population or number of locations from grid adjacency
#' @author Justin Moat. J.Moat@kew.org
#' @note Malin et al 2010 suggest cells width of 1/10 the longest distance. This uses the default grid (0,0)
#' @param thepoints dataframe of points of x,y
#' @param cellwidth in metres, (default = 1/10 the longest axis of the set of points)
#' @param neighborhood either rook (up,down,left and right) or queen (rook plus diagonals) default = queen
#' @param returnV three switches for return values 
#' S = simply the number of sub-pop or locations 
#' AREA = returns the area of the cells (use aoo for useful areas)
#' SF = return the simple feature for plotting and mapping
#' #' @return number of sub-pop/locations
#' @examples 
#'#Build some normally distributed point data around the Troodos mountains ~ 10 km diameter
#'thepoints <- normalofPs(50,0.1)
#'#shift to Troodos mountaions
#'thepoints <- data.frame(long = thepoints$X + 32.8794, lat = thepoints$Y + 34.9220)
#'#project the points
#'ppts <- simProjWiz(thepoints)
#'#get the number of sub-pop/locations using default method
#'subLocGrid(ppts)
#'#get area
#'subLocGrid(ppts,returnV="AREA")
#'#get the SF object and plot it
#'p <- subLocGrid(ppts,returnV="SF")
#'plot (p)
#'Using rook neighborhood
#'subLocGrid(ppts,neighborhood = "rook")
#'#with user defined cellwidth
#'subLocGrid(ppts,neighborhood = "rook",cellwidth = 2000)
#' @export
#' @import sf
#' @references 
#' Rivers, M. C., Bachman, S. P., Meagher, T. R., Lughadha, E. N., & Brummitt, N. A. (2010). Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data. Biodiversity and Conservation, 19(7), 2071-2085. https://link.springer.com/article/10.1007/s10531-010-9826-9
#' @seealso \code{\link{longestAxis}} to calculate the length of the longest axis
#' @seealso \code{\link{aoo}} to calculate simple AOO or
#' @seealso \code{\link{aooFixedGrid}} to calculate optimal AOO on a non-rotating grid or
#' @seealso \code{\link{aooFixedRotation}}to calculate optimal AOO with rotation and shift

subLocGrid <- function (thepoints,cellwidth=longestAxis(thepoints,returnV='S')/10,neighborhood="queen",returnV="S"){
  simpleaoopoly <- aoo(ppts,cellsize=cellwidth, returnV = "SF")
  if (neighborhood == "rook"){
    subp <- st_cast(st_union(simpleaoopoly),"POLYGON")
  } else {
    buf <- st_buffer(simpleaoopoly,cellwidth/10000)
    subp <- st_cast(st_union(buf),"POLYGON")
  }
  if(returnV=="SF"){return (subp)}
  if(returnV=="AREA"){
    parea <- sum(st_area(subp))
    units(parea) <- "km^2"
    return(parea)
  }
  length(subp)
}


##################################################################################
#calculates the number of Sub-population or number of locations from alpha hull method
##################################################################################
#' @title Sub-populations or locations using alpha hull method
#' @description 
#' Calculates the number of Sub-population or number of locations from alpha hull method
#' @author Justin Moat. J.Moat@kew.org
#' @note Malin et al 2010 suggest barrier width of 1/10 the longest distance. 
#' @param thepoints dataframe of points of x,y
#' @param barrierDis in metres, (default = 1/10 the longest axis of the set of points)
#' #' @param returnV four switches either 
#' S = for simply the number of sub-pop or locations 
#' AREA = returns the area of the cells (use aoo for useful areas)
#' SF = returns a multipolygon simple feature of the alpha hull, for mapping, plotting in ggplot or export to GIS systems
#' ALL = returns a multipolygon of all the triangular elements of the alpha hull
#' @return number of sub-pop/locations
#' @examples 
#'#Build some normally distributed point data around the Troodos mountains ~ 10 km diameter
#'thepoints <- normalofPs(50,0.1)
#'#shift to Troodos mountaions
#'thepoints <- data.frame(long = thepoints$X + 32.8794, lat = thepoints$Y + 34.9220)
#'#project the points
#'ppts <- simProjWiz(thepoints)
#'#get number of subpop with defaults
#'subLocAlpha(ppts)
#'#get area
#'subLocAlpha(ppts,returnV = "AREA")
#'#user defined distance
#'subLocAlpha(ppts,barrierDis=2000)
#'subpop <- subLocAlpha(ppts,returnV = "SF")
#'#plotting the results
#'library(ggplot2)
#'ggplot(data=subpop) + geom_sf() + geom_point(data=ppts,aes(X,Y))
#'
#' @export
#' @import sf
#' @references 
#' Rivers, M. C., Bachman, S. P., Meagher, T. R., Lughadha, E. N., & Brummitt, N. A. (2010). Subpopulations, locations and fragmentation: applying IUCN red list criteria to herbarium specimen data. Biodiversity and Conservation, 19(7), 2071-2085. https://link.springer.com/article/10.1007/s10531-010-9826-9
#' @seealso \code{\link{longestAxis}} to calculate the length of the longest axis
#' @seealso \code{\link{aHullMean}} to calculate the alpha hull using multiples of the mean lengths of the Delaunay triangulation
#' @note 
#' As the barrier distance increases for Alpha hulls, point become outliners with no area. This algorithm returns the number of groups and the number of these outliner for the sub population or location number
#' We would not recommend using Alpha hulls for sub populations or locations, but it is supplied here for consistency and for those that want to experiment.

subLocAlpha <- function (thepoints,barrierDis=longestAxis(thepoints,returnV='S')/10,returnV="S") {
  x <- cbind(thepoints$X,thepoints$Y)
  mp <- st_multipoint(x)
  mpm <- st_sfc(mp)
  st_crs(mpm)<- attr(thepoints,'crs')
  #build Delauney triangulation from sf package
  del <- st_triangulate(mp)
  #run through each triangle and get lengths
  t <- numeric()
  ind <- 0
  for (i in 1:length(del)){
    plist <- del[i][[1]][[1]]
    p1 <- st_point(plist[1,])
    p2 <- st_point(plist[2,])
    p3 <- st_point(plist[3,])
    t[3*i-2] <- st_distance(p1,p2)
    t[3*i-1] <- st_distance(p2,p3)
    t[3*i] <-  st_distance(p3,p1)
  }
  tridel_indexes <- unique(ceiling((which((t > barrierDis) %in% TRUE))/3))
  #if no triangle to remove then return normal convex hull
  if (length(tridel_indexes)==0) {
    ahullpoly <- st_union(del,del)
  } else {
    polys <- st_multipolygon(del[-tridel_indexes])
    ahullpoly <- st_union(polys,polys)
  }
  if (is.null(attr(thepoints,'crs'))){
    attr(thepoints,'crs') <- ''
    warning('Projection not set, it will be set to null')
  }
  ahullpoly <- st_sfc(ahullpoly)
  polys <- st_sfc(polys)
  st_crs(ahullpoly) <- attr(thepoints,'crs')
  st_crs(polys) <- attr(thepoints,'crs')
  #st_crs(ahullpoly, proj4text = attr(thepoints,'crs'))
  if (returnV=='SF'){return(ahullpoly)}
  if (returnV=='ALL'){return(polys)}
  if(returnV=="AREA"){
    parea <- st_area(ahullpoly)
    units(parea) <- "km^2"
    return(parea)
  }
  #get number of sub populations etc
  ip <- st_intersection (ahullpoly,mpm)
  #number of points used
  nopused <- length(st_geometry(ip)[[1]])/2
  #and therefore not used
  outliers <- nrow(thepoints) - nopused
  groups <- length(st_geometry(ahullpoly)[[1]])
  return(outliers + groups)
}


                 
