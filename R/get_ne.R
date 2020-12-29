#' @title Get Natural Earth country data
#'
#' @description Function to get the Natural Earth country boundary spatial data
#' @keywords Natural Earth
#' @import sf
#' @import rnaturalearth
#' @return sf object with data on global countries
#' @export

get_ne <- function() {
  countries10 <-
    ne_download(
      scale = 110,
      type = 'countries',
      category = 'cultural',
      returnclass = 'sf'
    )
  
}
