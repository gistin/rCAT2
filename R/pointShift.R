#' @title Point shift
#'
#' @description Shifts points by a set distance
#' @keywords shift
#' @return dataframe
#' @export

# Point shift - shift the point(s) randomly
pointShift <- function(points, distance){
  
  # generate a random point to buffer distance
  randPoint <- COfPs(1,distance)
  
  # add random distance to existing X and Y
  X = randPoint$x + points$X
  Y = randPoint$y + points$Y
  
  # make data frame
  shiftedPoints <- data.frame(X,Y)
}
