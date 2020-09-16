#utilities and extra scripts
#Thing about changing all the point tools to ptsSquare etc that way they are listed together


###########################################################
#calculates the longest axis from a set of points
###########################################################
#' @title Longest distance from a set of points
#' @description 
#' Calculates the longest distances from a set of points
#' @author Justin Moat. J.Moat@kew.org
#' @note Useful for scale for cellsize, Willis et al 2003 suggest 1/10 of this for cellsize for AOO calculations
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
#' @title Set of points in a square area
#' @description 
#' Builds a random set of points for a square area
#' @author Justin Moat. J.Moat@kew.org
#' @param nop number of points
#' @param gsize size of square (width or height), (default = 2000)
#' @return dataframe of points (x,y)
#' @examples
#' dfofpoints <- squareOfPs(100,1)
#' @export

squareOfPs <- function(noP,gsize=2000) {
  X <- runif (noP,0,gsize)
  Y <- runif (noP,0,gsize)
  df <- data.frame(X,Y)
  return (df)
}

###########################################################
#returns a set of random points for a oval/circle area    #
###########################################################
#' @title Set of a circle or oval points
#' @description 
#' Builds a random set of points for a circular or oval area
#' @author Justin Moat. J.Moat@kew.org
#' @param nop number of points
#' @param gsize size of square (width of longest size) (default 2000)
#' @param rot angle of rotation in radians (default 0.785398)
#' @param aspectRatio i.e. (major axis)/(minor axis), greater than 1, but if <1 it will just switch the axis. For a circle = 1
#' @return dataframe of points (x,y)
#' @examples
#' plot(ovalOfPs(100,1,1,0.5),asp=1)
#' @export

ovalOfPs <- function(noP,gsize=2000,rot=0.785398,aspectRatio=1) {
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
#' @title Set of annulus (doughnut) points
#' @description 
#' Builds random a set of random points for an annulus (doughnut) 
#' @author Justin Moat. J.Moat@kew.org
#' @param nop number of points
#' @param gsize scale of area (~ width of longest axis)
#' @param holes hole size 0 = none, 1 = ring of points (default=0.4)
#' @param aspectRatio i.e. (major axis)/(minor axis), greater than 1, but if <1 it will just switch the axis. For a circle = 1 (default=1)
#' @param rot angle of rotation in radians (default = 0)
#' @return dataframe of points (x,y)
#' @examples
#' doughnutofpoints <- doughnutOfPs (100,1,0.4,0.5,0.5)
#' plot(doughnutofpoints)
#' @export


doughnutOfPs <- function(noP=100,gsize=1,holes=0.4,aspectRatio=1,rot=0){
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
#returns a set of random points for normally distributed points
###########################################################
#' @title Set of normally distributed points
#' @description 
#' Builds random a set of normally distribute 
#' @author Justin Moat. J.Moat@kew.org
#' @param nop number of points
#' @param gsize scale of area (width of longest size), NB this is within 3 SD so some point may extend beyond this
#' @return dataframe of points (x,y)
#' @examples
#' normalpoints <- normalofPs(100,1)
#' plot(normalpoints)
#' @export

normalofPs <- function(nop=100,gsize=1){
  X<- rnorm(nop)*(gsize/3)
  Y<- rnorm(nop)*(gsize/3)
  df <- data.frame(X,Y)
  return (df)
}



