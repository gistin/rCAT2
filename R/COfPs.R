#' @title Circle of points
#'
#' @description Generate a circle of points
#' @keywords circle
#' @return dataframe
#' @export

# Circle of points - random point(s) to an error distance
COfPs <- function(noP,radius) {
  phi <- runif(noP,0,2*pi) #angle
  rho <- runif(noP,0,1) #distance
  x <- (sqrt(rho) * cos(phi)) * radius
  y <- (sqrt(rho) * sin(phi)) * radius
  df <- data.frame(x,y)
  return (df)
}
