###########################################################
#calculates the intial AOO, with simple grid 0,0 #
###########################################################
#' calculates a very simple AOO area from a set of points
#' @title Area of Occupancy (AOO)
#' @description 
#' Calculates the number of ocupied cells for Area of Occupancy from a set of points (x,y), usually in metres, with orgin 0,0
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints set of points in metres i.e. c(x,y)
#' @param cellsize size of cell (length) in metres
#' @return integer number of unique cells as an integer
#' @examples
#' x <- runif (20,0,10)
#' y <- runif (20,0,10)
#' df <- data.frame(x,y)
#' AOOsimp (df,2)
#' @seealso \code{\link{AOORating}} for AOO Ratings,
#' @export
#' @references
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Moat, J., Bachman, S., n.d. GeoCAT Geospatial Conservation Assessment Tool. URL http://geocat.kew.org/
AOOsimp <- function(thepoints,cellsize){
  bottomleftpoints <- floor(thepoints/cellsize)
  uniquecells <- unique(bottomleftpoints)
  #cellp <- data.frame(x=(uniquecells$x * cellsize), y=(uniquecells$y * cellsize))
  return(nrow(uniquecells))
}