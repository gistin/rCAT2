#data preparations scripts to help upload to IUCN etc

##########################################################
#Country list from a set of points in geographic coordinates#
##########################################################
#' calculates the MER of a set of numbers'
#' @title PRE-ALPHA Country list from a set of points in geographic coordinates
#' @description 
#' Calculates the country occurances from a set of points in lat long
#' @author Justin Moat. J.Moat@kew.org
#' @param thepoints dataframe of points ie c(x,y) in lat long
#' @return list of ISO3 codes
#' @examples
#' #built set of points
#' thepoints <- ptsNormal(29,0.1)
#'#shift to Peru, Colombia and Brazil tripoint
#' thepoints <- data.frame(long = thepoints$X - 69.9465972, lat = thepoints$Y - 4.2260167)
#' countrylist(thepoints)
#' @import sf
#' @import rnaturalearth
#' @export

countrylist <- function(thepointsll){
  #check to see if it's been downloaded
  if (!exists('countries10')){
    countries10 <- new.env()
    #note below does not always work, ie NE site can be down!!!
    countries10 <- ne_download(scale = 10, type='countries', category = 'cultural', returnclass='sf')
  }
  thepoints <- st_as_sf(thepoints, coords = c("long", "lat"), crs = 4326)
  clist <- st_intersects(thepoints,countries10)
  iso3l <- unique(counties10$ISO_A3[unlist(clist)])
  return(iso3l)
}
