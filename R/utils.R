#utilities and extra scripts
###########################################################
#calculates the longest axis from a set of points
###########################################################
#' @title Longest distance from a set of points
#' @description 
#' Calculates the longest distances from a set of points
#' @author Justin Moat. J.Moat@kew.org
#' @note Useful as a scale for cellsize and location buffers, Willis et al 2003 suggest 1/10 of this for cellsize for AOO calculations as does Rivers et al (2010) for buffer distance for sub-population or location calculations. 
#' @param thepoints dataframe of points of x,y
#' @param returnV, two switches either S for simply the distance or P for a dataframe of the two furthest points 
#' @return distance in metres or two points for the longest distance
#' @examples 
#' x <- runif (20,0,10)
#' y <- runif (20,0,10)
#' df <- data.frame(x,y) 
#' #distance only
#' longestAxis(df)
#' #two furthest points
#' dp <- longestAxis(df,'P')
#' plot(df, asp=1)
#' points(dp,pch=16)
#' @export
#' @references 
#' Willis, F., Moat, J., & Paton, A. (2003). Defining a role for herbarium data in Red List assessments: a case study of Plectranthus from eastern and southern tropical Africa. Biodiversity & Conservation, 12(7), 1537-1552.

longestAxis <- function (thepoints,returnV='S'){
  edgepoints <- thepoints[chull(thepoints),]
  distmax <- 0
  for (i in 1:(nrow(edgepoints)-1)){
    for (j in (i+1):nrow(edgepoints)){
      dist <- sqrt((edgepoints[i,1] - edgepoints[j,1])^2 + (edgepoints[i,2] - edgepoints[j,2])^2)
      if (dist > distmax){
        distmax <- dist
        p1 <- i
        p2 <- j
      }
    }
  }
  fpoints <- edgepoints[c(p1,p2),]
  if (returnV == "P"){return(fpoints)}
  distmax
}


##########################################################
#MER from a set of points#
##########################################################
#' calculates the MER of a set of numbers'
#' @title Minimum Enclosing Rectangle (mer)
#' @description 
#' Calculates the minimum enclosing rectangle (mer) from a set of points (x,y)
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points ie c(x,y)
#' @return vector of 4 doubles = xmin,xmax,ymin,ymax
#' @examples
#' x <- runif (20,0,10)
#' y <- runif (20,0,10)
#' df <- data.frame(x,y) 
#' mer (df)
#' @export



mer <- function(thepoints){
  xmin <- min(thepoints[1])
  xmax <- max(thepoints[1])
  ymin <- min(thepoints[2])
  ymax <- max(thepoints[2])
  return(c(xmin,xmax,ymin,ymax))
}

MER <- mer


###########################################################
#returns a set of random points for a square area         #
###########################################################
#' @title Builds a set of points in a square area
#' @description 
#' Builds a random set of points for a square area
#' @author Justin Moat. J.Moat@kew.org
#' @param nop number of points
#' @param gsize size of square (width or height), (default = 2000)
#' @return dataframe of points (x,y)
#' @examples
#' dfofpoints <- ptsSquare(100,1)
#' @export

ptsSquare <- function(noP,gsize=2000) {
  X <- runif (noP,0,gsize)
  Y <- runif (noP,0,gsize)
  df <- data.frame(X,Y)
  return (df)
}

###########################################################
#returns a set of random points for a oval/circle area    #
###########################################################
#' @title Builds a set of a circle or oval points
#' @description 
#' Builds a random set of points for a circular or oval area
#' @author Justin Moat. J.Moat@kew.org
#' @param nop number of points
#' @param gsize size of square (width of longest size) (default 2000)
#' @param rot angle of rotation in radians (default 0.785398)
#' @param aspectRatio i.e. (major axis)/(minor axis), can be greater than 1, but if so it will just switch the axis. (default =1 ie circle).
#' @return dataframe of points (x,y)
#' @examples
#' plot(ptsOval(100,1,1,0.5),asp=1)
#' @export

ptsOval <- function(noP,gsize=2000,rot=0.785398,aspectRatio=1) {
  xRatio <- sqrt((gsize ^ 2)/aspectRatio)
  yRatio <- xRatio * aspectRatio
  phi <- runif(noP,0,2*pi)
  rho <- runif(noP,0,1)
  x <- (sqrt(rho) * cos(phi)) * xRatio
  y <- (sqrt(rho) * sin(phi)) * yRatio
  df <- data.frame(x,y)
  rotationmatrix <- matrix(c(cos(rot),-sin(rot),sin(rot),cos(rot)), nrow = 2, ncol = 2)
  dfrot <- t(rotationmatrix %*% t(df))
  rdf <- as.data.frame(dfrot)
  colnames(rdf) <- c("X","Y")
  return (rdf)
}
###########################################################
#returns a set of random points for an annulus (doughnut)    #
###########################################################
#' @title Builds a set of annulus (doughnut) points
#' @description 
#' Builds random a set of random points for an annulus (doughnut shape) 
#' @author Justin Moat. J.Moat@kew.org
#' @param nop number of points
#' @param gsize scale of area (~ width of longest axis)
#' @param holes hole size, between 0 and 1, 0 = no hole, 1 = cirlce of points (default=0.4)
#' @param aspectRatio For an ellipse shape i.e. (major axis)/(minor axis), greater than 1, but if <1 it will just switch the axis. For a circle = 1 (default=1)
#' @param rot angle of rotation for the ellipse in radians (default = 0)
#' @return dataframe of points (x,y)
#' @examples
#' doughnutofpoints <- ptsDoughnut (100,1,0.4,0.5,0.5)
#' plot(doughnutofpoints)
#' @export


ptsDoughnut <- function(noP=100,gsize=1,holes=0.4,aspectRatio=1,rot=0){
  gsize <- 0
  theta<-runif(noP,0,2*pi)
  r<-sqrt(runif(noP,(holes*0.5)^2,0.5^2))
  ddf <-cbind((gsize+r*cos(theta))*aspectRatio,(gsize+r*sin(theta))/aspectRatio)
  rotationmatrix <- matrix(c(cos(rot),-sin(rot),sin(rot),cos(rot)), nrow = 2, ncol = 2)
  dfrot <- t(rotationmatrix %*% t(ddf))
  rdf <- as.data.frame(dfrot)
  colnames(rdf) <- c("X","Y")
  return(rdf)
}

###########################################################
#returns a set of random normally distributed points
###########################################################
#' @title Builds a set of normally distributed points
#' @description 
#' Builds random a set of normally distribute 
#' @author Justin Moat. J.Moat@kew.org
#' @param nop number of points
#' @param gsize scale of area (width of longest size), NB this is within 3 SD so some point may extend beyond this
#' @return dataframe of points (x,y)
#' @examples
#' normalpoints <- ptsNormal(100,1)
#' plot(normalpoints)
#' @export
ptsNormal <- function(nop=100,gsize=1){
  X<- rnorm(nop)*(gsize/3)
  Y<- rnorm(nop)*(gsize/3)
  df <- data.frame(X,Y)
  return (df)
}

#' Calculate Delauney triangles for alpha hull.
#' 
#' Performs Delauney triangulation on the provided points
#' and returns the triangles as a simple features collection.
#' 
#' @export
alphaTriangles <- function(x, y, crs) {
  
  if (is.null(crs)) {
    crs <- ""
  }
  
  #convert to matrix
  points <- cbind(x, y)
  points_geom <- st_multipoint(points)
  points_collection <- st_sfc(points_geom, crs=crs)
  
  if (is.na(st_crs(points_collection))) {
    warning("Missing CRS, assuming projection is in metres.")
  }
  
  #build Delauney triangulation from sf package
  del <- st_triangulate(points_collection)
  
  st_cast(st_sfc(del))
}

#' Calculate the distance between all vertices of a triangle
#' 
#' This is a utility to calculate the edge lengths of a triangle
#' in a Delauney triangulation.
#' 
#' @export
calculateDistances <- function(triangle) {
  points <- st_cast(st_sfc(triangle), "POINT")
  d_mat <- st_distance(points)
  d_mat[col(d_mat) == row(d_mat) + 1] 
}



