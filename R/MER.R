##########################################################
#MER from a set of points#
##########################################################
#' calculates the MER of a set of numbers'
#' @title Minimum Enclosing Rectangle (MER)
#' @description 
#' Calculates the minimum enclosing rectangle (MER) from a set of points (x,y)
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points ie c(x,y)
#' @return list of 4 doubles = xmin,xmax,ymin,ymax
#' @examples
#' x <- runif (20,0,10)
#' y <- runif (20,0,10)
#' df <- data.frame(x,y) 
#' MER (df)
#' @export



MER <- function(thepoints){
  xmin <- min(thepoints[1])
  xmax <- max(thepoints[1])
  ymin <- min(thepoints[2])
  ymax <- max(thepoints[2])
  return(c(xmin,xmax,ymin,ymax))
}