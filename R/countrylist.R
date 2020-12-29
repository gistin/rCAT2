#' Extract ISO country regions from a set of points in lat long
#' @param thepoints dataframe of points ie c(x,y) in lat long
#' @return dataframe or sf object with ISO country list
#' @examples
#' #built set of points
#' thepoints <- ptsNormal(29,0.1)
#'#shift to Peru, Colombia and Brazil tripoint
#' thepoints <- data.frame(long = thepoints$X - 69.9465972, lat = thepoints$Y - 4.2260167)
#' thepoints <- sf::st_as_sf(thepoints, coords = c("long", "lat"), crs = 4326)
#' countrylist(thepoints)
#' @import sf
#' @import dplyr
#' @export

countrylist <- function(thepointsll, returnV = 'S') {
  
  # read in the natural earth polygon data
  countries <- readRDS("data/countries_ne_10.rds")
  
  # intersect the points with the natural earth polygon polygons
  countries_points <- sf::st_intersection(thepointsll, countries)
  
  # df <- st_drop_geometry(df)
  countries_points_df <- sf::st_drop_geometry(countries_points)
  
  # remove duplicates from df
  countries_unique <- unique(countries_points_df)
  
  if (returnV == "SF") {
    
    # merge the table of unique with countries sf
    countries_unique_sf <- merge(countries, countries_unique)
    return(countries_unique_sf)
  }
  
  return(countries_unique)
}


