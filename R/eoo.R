#TODO
#sort out EOO min

##########################################################
#calculates a possible solution for min and max EOO from a set of points and error
##########################################################
#' @title BETA Minimum (near) Extent of Occurrence (EOO) from points with errors radius
#' @description 
#' Calculates the Extent of Occurrence in meters or returns a spatial polygon from a set of points (x,y)
#' Note this algorithm is in beta and will not always return the smallest area, this is true of very long and skinny sets of points. 
#' But will give a very quick approximation which in most situations  will suffice.
#' EOO's are constructed from the centroid and rays to each of the points on the convex hull and where these rays intercept the error circle
#' @note
#' BETA: it will go wrong with errors which overlap (i.e. small EOO)
#' Not quite min or max in most cases, but should be close.
#' Probably worth looking at when EOO is near thresholds.
#' @author Justin Moat. J.Moat@kew.org
#' @author Amelie Moat
#' @param thepoints dataframe of points in metres i.e. c(x,y,error or x,y)
#' @param errorfield field with the error radius if null, default error will be used
#' @param defaultRadius if no error field a default error will be applied to all points (default = 2000)
#' @param returnV switch to return different sets of results: \cr
#' S = Simple, returns just the minimum area in km2, (DEFAULT) \cr
#' EX = returns list for two areas -  near minimum and near maximum area for reference \cr
#' SF = returns a polygon simple feature for mapping, plotting in ggplot/plot or to export to a GIS format \cr 
#' SFA = returns a list of simple features with both the near minimum polygon and near maximum as a simple feature for mapping, plotting in ggplot or export to GIS systems
#' 
#' @return float_value area in km2 EOO polygon or sf polygon
#' @examples
#'#construct an oval of point for testing
#'thepoints <- ptsOval(19,0.5,deg2rad(45),0.5)
#'names(thepoints) <- c("lat","long")
#'#project the points
#'thepoints <- simProjWiz(thepoints)
#'#check crs
#'attr(thepoints,'crs')
#'#add some random errors between 0-5000m
#'thepoints$R <- runif(nrow(thepoints),0,5000)
#'#buffer them so you can view
#'ps <- st_as_sf(thepoints,coords=c('X','Y'))
#'psbuff <- st_buffer(ps,ps$R)
#'#normal EOO
#'eooPoly <- eoo(thepoints,"SF")
#'#set projection of bufferpoints
#'st_crs(psbuff) <- st_crs(eooPoly)
#'eooplot <- ggplot(data=eooPoly) + geom_sf(color="black",fill="green") + geom_sf(data=psbuff) 
#'eooplot
#'#normal eoo area
#'normaleookm <- eoo(thepoints)
#'#minimum eoo area
#'mineoo <- eooMin(thepoints,'R',,'S')
#'#percentage different
#'(normaleookm-mineoo)/normaleookm
#'# Get at sf object for plotting
#'mineoopoly <- eooMin(thepoints,'R',,"SF")
#'#add it to the gplot
#'eooplot <- eooplot + geom_sf(data=mineoopoly,fill=NA,col='red')
#'eooplot
#'#get the max as well
#'eooboth <- eooMin(thepoints,'R',,"SFA")
#'#add to plot
#'eooplot <- eooplot + geom_sf(data=eooboth$max,fill=NA,col='blue')
#'eooplot
#'###################
#'#with a content error radius
#'#construct an oval of point for testing
#'thepoints <- ptsOval(19,0.5,deg2rad(45),0.5)
#'names(thepoints) <- c("lat","long")
#'#project the points
#'thepoints <- simProjWiz(thepoints)
#'#
#'eooMin(thepoints,,2000,'S')
#'plot(eooMin(thepoints,,2000,'SF'))
#' @seealso \code{\link{ratingEoo}} for EOO Ratings
#' @export
#' @importFrom grDevices chull
#' @importFrom pracma polyarea
#' @import sf
#' @references
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591
#' 
eooMin <- function(thepoints, errorfield='R', defaultRadius=2000, returnV="S") {
  #check error field exists if not add default
  if (any(colnames(thepoints)==errorfield)){
    #rename the field for next bits
    colnames(thepoints)[which(colnames(thepoints)==errorfield)] <- "R"
  } else { #set default error
      thepoints$R <- defaultRadius
  }
  #get EOO as sf
  eooPoly <- eoo(thepoints,"SF")
  
  #get centroid as a matrix
  #NB need to check if centroid or true COG is better?
  midp <- as.numeric(st_coordinates(st_centroid(eooPoly)))
  #just get edge the points on the convex hull
  edgepoints <- thepoints[chull(thepoints),]
  #empty matrix to fill in loop
  innerPs <- matrix(ncol=2,nrow=nrow(edgepoints))
  outerPs <- matrix(ncol=2,nrow=nrow(edgepoints))
  #loop to get inner and outer intercepts
  for (i in 1:nrow(edgepoints)){
    ps <- l_c_intercepts(midp,c(edgepoints[i,]$X,edgepoints[i,]$Y),edgepoints[i,]$R)
    innerPs[i,] <- c(ps[1],ps[2])
    outerPs[i,] <- c(ps[3],ps[4])
  }
  # need to chull this set again in case a point can be dropped, as it may no longer be on the edge
  innerPs <- innerPs[chull(innerPs),]
  outerPs <- outerPs[chull(outerPs),]
  
  minArea <- -polyarea(innerPs[,1],innerPs[,2])/1000000
  maxArea <- -polyarea(outerPs[,1],outerPs[,2])/1000000
  
  if (returnV == "SF") {
    constructPolygon(innerPs[, 1], innerPs[, 2], attr(thepoints,'crs'))
  } else if (returnV == "EX") {
    list(
      min = minArea, 
      max = maxArea
    )
  } else if (returnV == "SFA"){
    minPoly <- constructPolygon(innerPs[, 1], innerPs[, 2], attr(thepoints,'crs'))
    maxPoly <- constructPolygon(outerPs[, 1], outerPs[, 2], attr(thepoints,'crs'))
    
    st_sf(
      measure=c("min", "max"),
      eoo=c(minArea, maxArea),
      geometry=c(minPoly, maxPoly)
    )
  } else {
    minArea
  }
}

##########################################################
#returns intercept for two points and radius
##########################################################
#example
#l_c_intercepts(c(2,1),c(-6,8),2)
l_c_intercepts <- function(midp,edgepoint,R){
  Ax <- midp[1]
  Ay <- midp[2]
  Bx <- edgepoint[1]
  By <- edgepoint[2]
  #compute the euclidean distance between A and B
  LAB <- sqrt((Bx-Ax)^2+(By-Ay)^2)
  #compute the direction vector D from A to B
  Dx <- (Bx-Ax)/LAB
  Dy = (By-Ay)/LAB
  # first intersection point
  Fx <- (LAB-R)*Dx + Ax
  Fy <- (LAB-R)*Dy + Ay
  # second intersection point
  Gx <- (LAB+R)*Dx + Ax
  Gy <- (LAB+R)*Dy + Ay
  return(c(Fx,Fy,Gx,Gy))
}

#' Construct a polygon from vertices.
#' 
#' Accepts an x and a y vector to define the vertices of
#' the polygon, to make it easier 
constructPolygon <- function(x, y, crs){
  points <- cbind(x, y)
  
  is_closed <- all(points[1,] == points[nrow(points),])
  
  if (! is_closed) {
    points <- rbind(points, points[1,])  
  }
  
  geom <- st_polygon(list(points))
  
  # put geometry into an sfc so we can attach a crs
  if (is.null(crs)) {
    crs <- ""
  }
  
  polygon <- st_sfc(geom, crs=crs)
  
  if (is.na(st_crs(polygon))) {
    warning("No valid CRS provided so setting it to `NA`")
  }
  
  polygon
}


###########################################################
#calculates the EOO area
###########################################################
#' @title Extent of Occurrence (EOO) Area
#' @description 
#' Calculates the Extent of Occurrence in km2 or returns a simple feature polygon from a set of points (x,y)
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points in metres i.e. c(x,y)
#' @param returnV switch to return different sets of results: \cr
#' S = Simple, returns just the minimum area in km2, (DEFAULT) \cr
#' SF = returns a polygon simple feature for mapping, plotting in ggplot/plot or to export to a GIS format
#' 
#' @return float_value area of EOO polygon or sf polygon
#' @note area returned is in x,y units, but negative as polygon is constructed anticlockwise
#' @examples
#' x <- runif (20,0,10)
#' y <- runif (20,0,10)
#' df <- data.frame(x,y) 
#' eoo (df)
#' #######
#' spoly <- eoo (df,TRUE)
#' plot(spoly)
#' @seealso \code{\link{ratingEoo}} for EOO Ratings
#' @export
#' @importFrom grDevices chull
#' @importFrom pracma polyarea
#' @import sf
#' @references
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591

eoo <- function(points, returnV="S") {
  if (! "X" %in% colnames(thepoints) | ! "Y" %in% colnames(thepoints)) {
    stop("Point coordinates must be supplied in columns named 'X' and 'Y'.")
  }
  
  hull_idx <- chull(points)
  hull <- points[hull_idx,]
  
  area <- polyarea(x=hull$X, y=hull$Y)
  # hull is constructed backwards, so area is negative and in m^2
  area <- -1 * area / 1e6
  
  if (returnV == "S") {
    area
  } else {
    constructPolygon(hull$X, hull$Y, attr(points, "crs"))
  }
} 

#depreciated but kept for compatabity with rCAT 1.6

EOOarea <- function(thepoints) {
  EOOpolyid <- chull(thepoints)
  EOOpp <- thepoints[EOOpolyid,]
  harea <- polyarea(x=EOOpp$x,y=EOOpp$y)
  #check for Area = NA ie when we only have one point
  if (is.na(harea)){harea <- 0}
  return(harea)
}

