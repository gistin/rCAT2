###########################################################
#Simple area Projection Wizard
###########################################################
#' Simple equal area projection wizard
#' @title Simple Projection Wizard
#' @description 
#' Projects any set of lat long points to a "suitable" area projection, based on it's "true centre of gravity"
#' @author Justin Moat. J.Moat@kew.org
#' @note
#' Based around a simple continental projection, using two sets of projections
#' equal area cylindrical = Cylindrical equal-area = 8287
#' equal area azimuthal for polar (above 70) = Lambert azimuthal equal-area
#'  
#' note these are not cartographically pleasing projections, they are just so we can get the data into something simple for areal analysis
#' See below for a more cartographically pleasing projection engine
#' 
#' Šavric, B., Jenny, B., Jenny, H., 2016. Projection Wizard – An Online Map Projection Selection Tool. Cartogr. J. 53, 1–9. doi:10.1080/00087041.2015.1131938
#' @param thepoints set of points in latitude and longtitude ie c(lat,long)
#' @param thecentre one point ie c(lat,long)
#' 
#' @return set of points in metres (x,y)
#' @examples 
#'lat <- runif (200,-24,-12)
#'long <- runif (200,43,51)
#'ll <- data.frame(lat,long)
#'cp <- trueCOGll(ll)
#'pointsprojected <- simProjWiz(ll,cp)
#' @references 
#' Šavric, B., Jenny, B., Jenny, H., 2016. Projection Wizard – An Online Map Projection Selection Tool. Cartogr. J. 53, 1–9. doi:10.1080/00087041.2015.1131938
#' 
#' Snyder, J.P., 1987. Map projections: A working manual, Professional Paper. Washington, D.C.
#' @export
#' @import sp
#' @import rgdal



######################################################################
#simple projection wizard#
######################################################################
#determining projection around on center of points
#based around a simplied continental scheme. Also see:
#Šavric, B., Jenny, B., Jenny, H., 2016. Projection Wizard – An Online Map Projection Selection Tool. Cartogr. J. 53, 
#1–9. doi:10.1080/00087041.2015.1131938
#two sets of projections
#equal area cylindrical = Cylindrical equal-area = 8287
#equal area azimuthal for polar (above 70) = Lambert azimuthal equal-area = 
#note these are not cartographically pleasing projections, they are just so we can get the data into something simple for areal analysis
######################################################################


simProjWiz <- function(thepoints,thecentre){
  #setup and set projection to WGS84
  coordinates(thepoints) <- c("long", "lat")
  proj4string(thepoints) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") # geographic and WGS 84
  #depending on centre point
  if((thecentre$lat < 70) & (thecentre$lat > -70)){
    CRSstring <- paste("+proj=cea +lon_0=", thecentre$long,   " +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep = "")
  } else {
    CRSstring <- paste("+proj=laea +lat_0=", thecentre$lat," +lon_0=", thecentre$long, " +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",sep = "")
  }
  CRS.new <- CRS(CRSstring)
  #reproject
  xysp <- spTransform(thepoints, CRS.new)
  xy <- as.data.frame(xysp)
  #rename to x and y as not longer lat long
  colnames (xy) <- c("x","y")
  return(xy)
}

######################################################################
#calculates 'true' centre of gravity from a set of lat long points in decimal degrees   #
#note z from mean of cartesian give some idea of weighted spread on the globe#
######################################################################
#' @title True centre of gravity from a set of Lat longs
#' @description 
#' Calculates the "true" centre of gravity (weighted) from a set of lat longs, using cartesian geometry
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints set of points c(lat,long)
#' @return a point (lat,long) from centre
#' @examples 
#'lat <- runif (200,-24,-12)
#'long <- runif (200,43,51)
#'ll <- data.frame(lat,long)
#'cp <- trueCOGll(ll)
#' @references Descartes, R., 1637. Discours de la methode. A Leyde, De l’imprimerie de I. Maire, Paris.
#' @export


trueCOGll <-function(thepoints){
  
  llrad <- deg2rad(thepoints) #to radians
  cartp <- ll2cart(llrad$lat,llrad$long) #to cartesian
  mp <- data.frame(x=mean(cartp$x),y=mean(cartp$y),z=mean(cartp$z)) #central point
  pmp <- pro2sph(mp$x,mp$y,mp$z) #projection to surface
  pmprll <- cart2ll(pmp$x,pmp$y,pmp$z) #to ll in radians
  pmpll <- rad2deg(pmprll) #to degrees
  return(data.frame(lat=pmpll$latr,long=pmpll$longr))
  
}







######################################################################
#calculates the Cartesian cordinates (x,y,z) from lat long in radians#
######################################################################
#' @title Geographic coordinates to cartesian (x,y,z)
#' @description 
#' Calculates the Cartesian cordinates (x,y,z) from lat long in radians
#' @author Justin Moat. J.Moat@kew.org
#' @param latr latitude point in radians
#' @param longr longtitude point in radians
#' @return dataframe of x,y,z
#' @examples 
#'lat <- runif (200,-24,-12)
#'long <- runif (200,43,51)
#'thepoints <- data.frame(lat,long)
#'llrad <- deg2rad(thepoints)
#'cartp <- ll2cart(llrad$lat,llrad$long)
#' @references Descartes, R., 1637. Discours de la methode. A Leyde, De l’imprimerie de I. Maire, Paris.
#' @export

ll2cart <- function(latr,longr){
  x <- cos(latr) * cos(longr)
  y <- cos(latr) * sin(longr)
  z <- sin(latr)
  return(data.frame(x,y,z))
}

######################################################################
#calculates the lat long cordinates in radians from Cartesian (x,y,z)#
######################################################################
#' @title Cartesian (x,y,z) to Geographic coordinates
#' @description 
#' calculates the latitude and longtitude cordinates in radians from Cartesian coordinates (x,y,z)
#' @author Justin Moat. J.Moat@kew.org
#' @param x East to West coordinate in metres
#' @param y South to North coordinate in metres
#' @param z height coordinate in metres
#' @return dataframe of latitude,longtitude
#' @export

cart2ll <-function (x,y,z){
  latr <- asin(z)
  longr <- atan2(y,x)
  return(data.frame(latr,longr))
}


######################################################################
#calculates Cartesian (x,y,z), projected from the centre of the sphere 
#to the earth surface, returns cartesian (x,y,z)
#used to calculate "true" centre of set of lat longs
# http://stackoverflow.com/questions/9604132/how-to-project-a-point-on-to-a-sphere
######################################################################
#' @title Cartesian coordinate projection
#' @description 
#' calculates Cartesian (x,y,z), projected from the centre of the sphere 
#' to the earth surface, returns cartesian (x,y,z)
#' used to calculate "true" centre of set of lat longs
#' @author Justin Moat. J.Moat@kew.org
#' @note
#' http://stackoverflow.com/questions/9604132/how-to-project-a-point-on-to-a-sphere
#' @param x East to West coordinate in metres
#' @param y South to North coordinate in metres
#' @param z height coordinate in metres
#' @return x,y,z
#' @references Descartes, R., 1637. Discours de la methode. A Leyde, De l’imprimerie de I. Maire, Paris.
#' @export

pro2sph <- function (x,y,z){
  sc <- 1/sqrt(x^2 + y^2 + z^2)
  x <- x * sc
  y <- y * sc
  z <- z * sc
  return(data.frame(x,y,z))
}

######################################################################
#radians to degrees and degrees to radians
######################################################################
#' @title Radians to Degrees
#' @description 
#' Calculates radians from degrees or degrees from radians
#' @author Justin Moat. J.Moat@kew.org
#' @param rad number in radians
#' @return number
#' @examples 
#' b <- 0.392699
#' rad2deg(b)
#' @export

rad2deg <- function(rad) {(rad * 180) / (pi)}

######################################################################
#radians to degrees and degrees to radians
######################################################################
#' @title 
#' Degrees to radians
#' @description 
#' Calculates radians from degrees or degrees from radians
#' @author Justin Moat. J.Moat@kew.org
#' @param deg number in degrees
#' @return number
#' @examples 
#' a <- 30
#' deg2rad(a)
#' @export

deg2rad <- function(deg) {(deg * pi) / (180)}