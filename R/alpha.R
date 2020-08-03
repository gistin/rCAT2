###########################################################
#calculates the alpha hull from a set of points
###########################################################
#' @title Alpha hull 
#' @description
#' Calculates the Alpha hull from a set of points and m (multiple of the mean distance between points)
#' Processing time for large dataset will be slow. uses sf's st_triangulate initially to build Delauney triangulation of all points 
#' Algorithms derived from and for details of use see: 
#' Burgman, M. A., & Fox, J. C. (2003). Bias in species range estimates from minimum convex polygons: implications for conservation and options for improved planning. Animal Conservation, 6(1), 19-28.
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points in metres i.e. c(x,y)
#' @param multiple multiple of the mean distance between point pairs, used to remove triangles from the convex hull (default=2)
#' @param returnV, switches to return different sets of results: 
#' S = Simple, returns just the area in km2, (DEFAULT)
#' SF = returns a multipolygon simple feature of the alpha hull, for mapping, plotting in ggplot or export to GIS systems
#' ALL = returns a multipolygon of all the triangular elements of the alpha hull
#' @return dependent on switch, default is area in km2
#' @note 
#' @examples
#'#Build and project some points
#'thepoints <- squareOfPs(19,0.1)
#'names(thepoints) <- c("lat","long")
#'thepoints <- simProjWiz(thepoints)
#'attr(thepoints,'crs')
#'
#'#just return area in km2
#'aHullMean(thepoints,2)
#plot the alpha hull
#'plot(aHullMean(thepoints,2,'SF'))
#'points(thepoints)
#plot the Delauney triangulation alpha hull
#'plot(aHullMean(thepoints,2,'SFE'),col='grey')
#'points(thepoints)
#' @export
#' @import sf
#' @references
#' Burgman, M. A., & Fox, J. C. (2003). Bias in species range estimates from minimum convex polygons: implications for conservation and options for improved planning. Animal Conservation, 6(1), 19-28.



aHullMean <- function(thepoints,multiple=2,returnV='S') {
  #convert to matrix
  x <- cbind(thepoints$X,thepoints$Y)
  mp <- st_multipoint(x)
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
  #calculate the mean
  meand <- mean(t)
  #work out the triangles that are going to be deleted
  #get the index of the polygons which are greater
  #NB /3 and ceiling to get to index of triangles
  tridel_indexes <- unique(ceiling((which((t > meand*multiple) %in% TRUE))/3))
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
  return(as.numeric(st_area(ahullpoly)/1000000))
}
