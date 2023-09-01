#TD large number of multiple  returns odd results
#TD mean calculation uses repeated internal triangle lengths, does this matter???? ie external lengths only occur once internal twice
###########################################################
#calculates the alpha hull from a set of points
###########################################################
#' @title Alpha hull 
#' @description
#' Calculates the Alpha hull from a set of points and (multiple of the mean distance between points)
#' Processing time for large dataset will be slow. Uses sf's st_triangulate initially to build Delaunay triangulation of all points. 
#' Algorithms derived from and for details of use see: 
#' Burgman, M. A., & Fox, J. C. (2003). Bias in species range estimates from minimum convex polygons: implications for conservation and options for improved planning. Animal Conservation, 6(1), 19-28.
#' @author Justin Moat. J.Moat@kew.org Steve Bachman S.Bachman@kew.org
#' @param thepoints dataframe of points in metres i.e. c(x,y)
#' @param multiple multiple of the mean distance between point pairs, used to remove triangles from the set of Delaunay triangulations (default=2)
#' @param returnV, switches to return different sets of results: \cr
#' S = Simple, returns just the area in km2, (DEFAULT) \cr
#' SF = returns a multipolygon simple feature of the alpha hull, for mapping, plotting in ggplot or export to GIS systems \cr
#' ALL = returns a multipolygon of all the triangular elements of the alpha hull \cr
#' @return dependent on switch, default is area in km2 
#' @note none
#' @examples
#'#Build and project some points
#'thepoints <- ptsSquare(19,0.1)
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
#'
aHullMean <- function(thepoints, multiple=2, returnV='S') {
  crs <- attr(thepoints, "crs")
  
  triangles <- alphaTriangles(thepoints$X, thepoints$Y, crs=crs)
  
  distances <- lapply(triangles, calculateDistances)
  
  #calculate the overall mean
  mean_distance <- mean(do.call(c, distances))
  threshold <- mean_distance * multiple
  
  # get logical mask of triangles to keep
  keep_triangle <- sapply(distances, function(x) all(x <= threshold))
  
  if (returnV == 'SF') {
    st_union(triangles[keep_triangle])
  } else if (returnV == 'ALL') {
    polys <- st_multipolygon(triangles[keep_triangle])
    st_sfc(polys, crs=crs)
  } else {
    # assumes CRS is set
    hull <- st_union(triangles[keep_triangle])
    area <- st_area(hull) 
    
    if (is.null(attr(area, "units"))) {
      area <- area / 1e6
    }
    
    units(area) <- "km^2"
    area
  }
}

