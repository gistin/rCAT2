##########################################################
#calculates a possible solution for min and max EOO from a set of points and error
##########################################################
#' @title Minimum (near) Extent of Occurrence (EOO) from points with errors
#' @description 
#' Calculates the Extent of Occurrence in meters or returns a spatial polygon from a set of points (x,y)
#' Note this algorithm will not always return the smallest area, this is true of very long and skinny sets of points. 
#' But will give a very quick approximation which in most situations  will work.
#' EOO's are constructed from the centroid and rays to each of the points on the convex hull and where these rays intercept the error circle
#' @note
#' BETA: it will go wrong with errors which overlap (ie small EOO)
#' Not quite min or max in most cases, but close
#' 
#' @author Justin Moat. J.Moat@kew.org
#' @author Amelie Moat
#' @param thepoints dataframe of points in metres i.e. c(x,y,error or x,y)
#' @param errorfield field with the error radius if null, default error will be used
#' @param defaultRadius if no error field a default error will be applied to all points (default = 2000)
#' @param returnV switch to return different sets of results: 
#' S = Simple, returns just the minimum area in km2, (DEFAULT)
#' EX = returns list for two areas -  near minimum and near maximum area for reference
#' SF = returns the minimum polygon as a simple feature for mapping, plotting in ggplot or export to GIS systems
#' SFA = returns a list with both the near minimum polygon and near maximum as a simple feature for mapping, plotting in ggplot or export to GIS systems
#' 
#' @return float_value area in km2 EOO polygon or sf polygon
#' @examples
#'#construct an oval of point for testing
#'thepoints <- ovalOfPs(19,0.5,deg2rad(45),0.5)
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
#'thepoints <- ovalOfPs(19,0.5,deg2rad(45),0.5)
#'names(thepoints) <- c("lat","long")
#'#project the points
#'thepoints <- simProjWiz(thepoints)
#'#
#'eooMin(thepoints,,2000,'S')
#'plot(eooMin(thepoints,,2000,'SF'))
#' @seealso \code{\link{eooRating}} for EOO Ratings
#' @export
#' @importFrom grDevices chull
#' @importFrom pracma polyarea
#' @import sf
#' @references
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591
#' 
eooMin <- function(thepoints,errorfield='R',defaultRadius=2000,returnV="S") {
  #check error field exists if not add default
  if (any(colnames(thepoints)==errorfield)){
    #rename the field for next bits
    colnames(thepoints)[which(colnames(thepoints)==errorfield)] <- "R"
  } else { #set default error
      thepoints$R <- defaultRadius
  }
  #get EOO as sf
  eooPoly <- eoo(thepoints,"SF")
  #just stops a warning, not sure why
  st_agr(eooPoly) = "constant"
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
  ######need to chull this set again in case a point can be dropped, as it may no longer be on the edge
  innerPs <- innerPs[chull(innerPs),]
  outerPs <- outerPs[chull(outerPs),]
  ###

  
  if (returnV == "SF") {return(polyCon(innerPs,attr(thepoints,'crs')))}
  if (returnV == "EX"){return(list(min = -polyarea(innerPs[,1],innerPs[,2])/1000000, max = -polyarea(outerPs[,1],outerPs[,2])/1000000))}
  if (returnV == "SFA"){
    minpoly <- polyCon(innerPs,attr(thepoints,'crs'))
    maxpoly <- polyCon(outerPs,attr(thepoints,'crs'))
    return(list(min=minpoly,max=maxpoly))}
  -polyarea(innerPs[,1],innerPs[,2])/1000000
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

##########################################################
#eoo polygon constructor
#########################################################
polyCon <- function(pointsmatrix,crs){
  pointsdf <- data.frame(X=pointsmatrix[,1],Y=pointsmatrix[,2],id=1)
  xys <- st_as_sf(pointsdf,coords=c("X","Y"))
  poly <- st_sf(
    aggregate(
      xys$geometry,
      list(xys$id),
      function(g){
        st_cast(st_combine(g),"POLYGON")
      }
    ))
  #set crs
  st_crs(poly) <- crs
  return(poly)
}


###########################################################
#calculates the EOO area
###########################################################
#' @title Extent of Occurrence (EOO) Area
#' @description 
#' Calculates the Extent of Occurrence in km2 or returns a spatial polygon from a set of points (x,y)
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points in metres i.e. c(x,y)
#' @param returnV switch to return different sets of results: 
#' S = Simple, returns just the minimum area in km2, (DEFAULT)
#' SF = returns a polygon simple feature for mapping, plotting in ggplot or export to GIS systems
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
#' @seealso \code{\link{eooRating}} for EOO Ratings
#' @export
#' @importFrom grDevices chull
#' @importFrom pracma polyarea
#' @import sf
#' @references
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591

eoo <- function(thepoints,returnV="S") {
  EOOpolyid <- chull(thepoints)
  EOOpp <- thepoints[EOOpolyid,]
  EOOpp$id <- 1
  #construct polygons
  if (returnV=="S"){return(-polyarea(x=EOOpp$X,y=EOOpp$Y)/1000000)
    } else {
    xys <- st_as_sf(EOOpp,coords=c("X","Y"))
    poly <- st_sf(
      aggregate(
        xys$geometry,
        list(xys$id),
        function(g){
          st_cast(st_combine(g),"POLYGON")
        }
      ))
    if (is.null(attr(thepoints,'crs'))){
      attr(thepoints,'crs') <- ''
      warning('Projection not set, it will be set to null')
    }
    st_crs(poly) <- attr(thepoints,'crs')
    return(poly)
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

