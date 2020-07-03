###########################################################
#calculates the mode of a set of data - note the capital M#
###########################################################
#from http://stackoverflow.com/questions/2547402/standard-library-function-in-r-for-finding-the-mode
#' calculates the mode of a set of numbers'
#' @title Mode of a test of points
#' @description 
#' Calculates the mode (the value that occurs most often) of a set of data - note the capital M
#' @author Justin Moat. J.Moat@kew.org
#' @note
#' Orginally from http://stackoverflow.com/questions/2547402/standard-library-function-in-r-for-finding-the-mode 
#' @param x set of numbers
#' @return number
#' @examples 
#' a <- c(5,5,5,6,7,8,9)
#' Mode(a)
#' @export

Mode <- function(x) { 
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}