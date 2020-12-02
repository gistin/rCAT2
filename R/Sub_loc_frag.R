#Scripts for Sub-population and number of locations (same methods) or fragmentation
##################################################################################
#calculates the number of Sub-population or number of locations from buffer method
##################################################################################
#' @title Sub-population or number of locations using the buffer method
#' @description 
#' Calculates the number of Sub-population or number of locations using the buffer method
#' @author Justin Moat. J.Moat@kew.org
#' @note Malin et al 2010 suggest buffer size of 1/10 the longest distance
#' @param thepoints dataframe of points of x,y
#' @param bufferradius in metres, (default = 1/10 the longest axis of the set of points)
#' @param returnV three switches for return values \cr
#' S = simply the number of sub-pop or locations \cr
#' AREA = returns the area of the buffered output \cr
#' SF = return the simple feature for plotting, mapping and export to a GIS
#' @return number of sub-pop/locations
#' @examples 
#'#Build some normally distributed point data around the Troodos mountains ~ 10 km diameter
#'thepoints <- ptsNormal(50,0.1)
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

subLocBuf <- function (thepoints, bufferradius=longestAxis(thepoints,returnV='S')/10, returnV="S"){
  thepointsSF <- st_as_sf(thepoints, coords = c("X", "Y"), crs=attr(thepoints, 'crs'))
  buf <- st_union(st_buffer(thepointsSF, bufferradius))
  
  if (returnV == "AREA") {
    area <- st_area(buf)
    
    units(area) <- "km^2"
    area
  } else if (returnV == "SF") {
    st_cast(buf, "POLYGON")
  } else {
    length(st_cast(buf, "POLYGON"))
  }
}

##################################################################################
#calculates the number of Sub-population or number of locations from cell adjacency method
##################################################################################
#' @title Sub-population or number of locations using grid adjacency method
#' @description 
#' Calculates the number of Sub-population or number of locations using the grid adjacency method
#' @author Justin Moat. J.Moat@kew.org
#' @note Malin et al 2010 suggest cells width of 1/10 the longest distance. This uses the default grid (0,0)
#' @param thepoints dataframe of points of x,y
#' @param cellwidth in metres, (default = 1/10 the longest axis of the set of points)
#' @param neighborhood either rook (up,down,left and right) or queen (rook plus diagonals) default = queen
#' @param returnV three switches for return values \cr
#' S = simply the number of sub-pop or locations \cr
#' AREA = returns the area of the cells (use AOO for useful areas methods for IUCN assessments) \cr
#' SF = return the simple feature for plotting, mapping and export to a GIS
#' #' @return number of sub-pop/locations
#' @examples 
#'#Build some normally distributed point data around the Troodos mountains ~ 10 km diameter
#'thepoints <- ptsNormal(50,0.1)
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

subLocGrid <- function (thepoints, cellwidth=longestAxis(thepoints,returnV='S')/10, neighborhood="queen", returnV="S"){
  simpleaoopoly <- aoo(thepoints, cellsize=cellwidth, returnV = "SF")
  if (neighborhood == "rook"){
    subp <- st_union(simpleaoopoly)
  } else {
    buf <- st_buffer(simpleaoopoly, cellwidth/10000)
    subp <- st_union(buf)
  }
  
  if (returnV == "AREA") {
    area <- st_area(subp)
    units(area) <- "km^2"
    area
  } else if (returnV == "SF") {
    st_cast(subp, "POLYGON")
  } else {
    length(st_cast(subp, "POLYGON"))
  }
}


##################################################################################
#calculates the number of Sub-population or number of locations from alpha hull method
##################################################################################
#' @title Sub-populations or locations using the alpha hull method
#' @description 
#' Calculates the number of Sub-population or  locations using the alpha hull method
#' @author Justin Moat. J.Moat@kew.org
#' @note Malin et al 2010 suggest barrier width of 1/10 the longest distance. 
#' @param thepoints dataframe of points of x,y
#' @param barrierDis in metres, (default = 1/10 the longest axis of the set of points)
#' @param returnV four switches either \cr
#' S = for simply the number of sub-pop or locations \cr
#' AREA = returns the area of the cells (use AOO or EOO for useful areas for IUCN assessment) \cr
#' SF = returns a multipolygon simple feature of the alpha hull, for mapping, plotting in ggplot or export to a GIS  \cr
#' ALL = returns a multipolygon of all the triangular elements of the alpha hull
#' @return number of sub-pop/locations
#' @examples 
#'#Build some normally distributed point data around the Troodos mountains ~ 10 km diameter
#'thepoints <- ptsNormal(50,0.1)
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
#' As the barrier distance increases for Alpha hulls, point become outliers with no area. This algorithm returns the number of grouped points (in the hull) and the number of these outlying points for the sub population or location calculations. 
#' We would not recommend using Alpha hulls for sub populations or locations, but it is supplied here for consistency and for those that want to experiment.

subLocAlpha <- function (thepoints, barrierDis=longestAxis(thepoints,returnV='S')/10, returnV="S") {
  crs <- attr(thepoints, "crs")
  
  triangles <- alphaTriangles(thepoints$X, thepoints$Y, crs=crs)
  
  distances <- lapply(triangles, calculateDistances)
  
  # get logical mask of triangles to keep
  keep_triangle <- sapply(distances, function(x) all(x <= barrierDis))
  
  if (returnV == 'SF') {
    st_union(triangles[keep_triangle])
  } else if (returnV == 'ALL') {
    polys <- st_multipolygon(triangles[keep_triangle])
    st_sfc(polys, crs=crs)
  } else if (returnV == "AREA") {
    hull <- st_union(triangles[keep_triangle])
    area <- st_area(hull) 
    
    if (is.null(attr(area, "units"))) {
      area <- area / 1e6
    }
    
    units(area) <- "km^2"
    area
  } else {
    # number of subpopulations is the number of polygons
    # plus number of points outside of hull
    
    hull <- st_union(triangles[keep_triangle])
    point_sf <- st_as_sf(thepoints, coords=c("X", "Y"), crs=crs)
    
    in_hull <- st_intersects(point_sf, hull, sparse=FALSE)
    
    outliers <- sum(! in_hull)
    groups <- length(hull)
    
    outliers + groups
  }
}


                 
