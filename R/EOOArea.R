###########################################################
#calculates the EOO area
###########################################################
#' calculates the EOO area of a set of popints'
#' @title Extent of Occurance (EOO) Area
#' @description 
#' Calculates the Extent of Occurance from a set of points (x,y)
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points i.e. c(x,y)
#' @return float_value area of EOO polygon
#' @note area returned is in x,y units, but negative as polygon is constructed anticlockwise
#' @examples
#' x <- runif (20,0,10)
#' y <- runif (20,0,10)
#' df <- data.frame(x,y) 
#' EOOarea (df)
#' @seealso \code{\link{EOORating}} for EOO Ratings
#' @export
#' @importFrom grDevices chull
#' @importFrom pracma polyarea
#' @references
#' Bachman, S., Moat, J., Hill, A.W., de Torre, J., Scott, B., 2011. Supporting Red List threat assessments with GeoCAT: geospatial conservation assessment tool. Zookeys 126, 117–26. doi:10.3897/zookeys.150.2109 
#' 
#' Moat, J., Bachman, S., n.d. GeoCAT Geospatial Conservation Assessment Tool [WWW Document]. URL http://geocat.kew.org/
#' 
#' Joppa, L.N., Butchart, S.H.M., Hoffmann, M., Bachman, S.P., Akçakaya, H.R., Moat, J.F., Böhm, M., Holland, R.A., Newton, A., Polidoro, B., Hughes, A., 2016. Impact of alternative metrics on estimates of extent of occurrence for extinction risk assessment. Conserv. Biol. 30, 362–370. doi:10.1111/cobi.12591


EOOarea <- function(thepoints) {
  EOOpolyid <- chull(thepoints)
  EOOpp <- thepoints[EOOpolyid,]
  harea <- polyarea(x=EOOpp$x,y=EOOpp$y)
  #check for Area = NA ie when we only have one point
  if (is.na(harea)){harea <- 0}
  return(harea)
}